#!/usr/bin/env python3

import argparse
import logging
import re
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from google.cloud import storage
from google.cloud.storage import transfer_manager

logger = logging.getLogger("build_and_prep")

def setup_logging(log_file: Path) -> None:
    """Logs to stdout and to log_file"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode="w")],
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Builds, prepares, and validates submission files for seqsender")
    parser.add_argument("--metadata_file", type=Path, required=True, help="Terra table TSV, named <table_name>.tsv")
    parser.add_argument("--config_file", type=Path, required=True, help="Path to the seqsender configuration file")
    parser.add_argument("--sample_ids", nargs="+", required=True, help="Sample IDs to be submitted")
    parser.add_argument("--bioproject", type=str, required=True, help="NCBI BioProject ID to submit to")

    # Terra lists these in the same order as the sample IDs
    parser.add_argument("--assembly_fastas", nargs="*", type=Path, default=[], help="Required for GenBank; one assembly FASTA per sample")
    parser.add_argument("--read1_fastqs", nargs="*", type=Path, default=[], help="Required for SRA; one 'read1' FASTQ per sample")
    parser.add_argument("--read2_fastqs", nargs="*", type=Path, default=[], help="Optional for SRA; leave out for single-end reads; one 'read2' FASTQ per sample")

    # with a transfer bucket, NCBI fetches the reads from the bucket ('cloud'); without one, seqsender uploads them over FTP ('local')
    parser.add_argument("--sra_transfer_bucket", type=str, help="Optional; NCBI-readable bucket the passing samples' reads are copied to")

    parser.add_argument("--organism", type=str, required=True, help="seqsender organism; options include 'FLU', 'COV', or 'OTHER'")
    parser.add_argument("--biosample", action="store_true")
    parser.add_argument("--sra", action="store_true")
    parser.add_argument("--genbank", action="store_true")
    parser.add_argument("--validation_only", action="store_true", help="Validate only; don't copy the reads to --sra_transfer_bucket")
    parser.add_argument("--output_dir", type=Path, default=Path.cwd(), help="Directory to write the submission files to; populates 'builds' and 'submit' directories")
    return parser.parse_args()


@dataclass
class Sample:
    """A sample's staged files and its seqsender prep result"""
    sample_id: str
    sample_dir: Path  # builds/<sample_id>
    fasta: Optional[Path] = None  # builds/<sample_id>/<sample_id>.assembly.fasta
    reads: list[Path] = field(default_factory=list)  # in builds/raw_reads, read 1 first
    metadata: Optional[pd.DataFrame] = None  # the sample's Terra row plus the columns build_and_prep fills itself
    passed: bool = False
    reason: str = ""

    @property
    def metadata_csv(self) -> Path:
        return self.sample_dir / f"{self.sample_id}.metadata.csv"

    @property
    def prep_log(self) -> Path:
        return self.sample_dir / f"{self.sample_id}.seqsender_prep.log"


def check_inputs(args) -> list[str]:
    """
    Exits on problems with the inputs as a whole; problems with a single sample only fail that sample.
    Returns the seqsender database flags, e.g. ["--biosample", "--sra"]
    """
    db_flags = [f"--{name}" for name in ("biosample", "sra", "genbank") if getattr(args, name)]
    if not db_flags:
        raise SystemExit("ERROR: choose at least one of --biosample, --sra or --genbank")
    if args.sra and not args.read1_fastqs:
        raise SystemExit("ERROR: --read1_fastqs is required for SRA submission")
    if args.genbank and not args.assembly_fastas:
        raise SystemExit("ERROR: --assembly_fastas is required for GenBank submission")
    if args.genbank and args.organism not in ("FLU", "COV"):
        raise SystemExit(f"ERROR: GenBank submission is only supported for FLU and COV, not {args.organism}")
    if args.sra_transfer_bucket:
        if not args.sra_transfer_bucket.startswith("gs://"):
            raise SystemExit(f"ERROR: --sra_transfer_bucket must be a gs:// bucket, not {args.sra_transfer_bucket}")
        args.sra_transfer_bucket = args.sra_transfer_bucket.rstrip("/")

    for name in ("assembly_fastas", "read1_fastqs", "read2_fastqs"):
        files = getattr(args, name)
        if files and len(files) != len(args.sample_ids):
            raise SystemExit(f"ERROR: {len(files)} {name} given for {len(args.sample_ids)} sample_ids; they must be the same length and order")

    # read1/2 FASTQ files must be unique so they don't get overwritten in 'builds/raw_reads'
    if args.sra:
        names = [read.name for read in args.read1_fastqs + args.read2_fastqs]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise SystemExit(f"ERROR: read files must have unique names; repeated: {', '.join(duplicates)}")
    return db_flags


def stage_samples(
    builds_dir: Path,
    sample_ids: list[str],
    assembly_fastas: list[Path],
    read1_fastqs: list[Path],
    read2_fastqs: list[Path],
) -> list[Sample]:
    """
    Creates a 'builds/<sample>' directory per `sample_id` and moves:
      - the corresponding `assembly_fasta` input into 'builds/<sample>/<sample>.assembly.fasta'
      - the corresponding `read1/2_fastqs` input into 'builds/raw_reads/' (required for `seqsender prep`)

    Terra lists the corresponding input files in the same order as the sample IDs.

    Returns one `Sample` per `sample_id`.
    """
    raw_reads_dir = builds_dir / "raw_reads"
    raw_reads_dir.mkdir(parents=True, exist_ok=True)

    all_samples = []
    for i, sample_id in enumerate(sample_ids):
        sample = Sample(
            sample_id=sample_id,
            sample_dir=builds_dir / sample_id
        )
        # creating 'builds/<sample>' directory
        sample.sample_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created sample directory: {sample.sample_dir}")

        # inputs are moved rather than copied so the FASTQs aren't stored twice on the task's disk
        if assembly_fastas:
            sample.fasta = sample.sample_dir / f"{sample_id}.assembly.fasta"
            logger.info(f"Moving {assembly_fastas[i]} to {sample.fasta}")
            shutil.move(assembly_fastas[i], sample.fasta)

        # FASTQ read files keep their original file names
        for read_fastqs in (read1_fastqs, read2_fastqs):
            if read_fastqs:
                read = raw_reads_dir / read_fastqs[i].name
                logger.info(f"Moving {read_fastqs[i]} to {read}")
                shutil.move(read_fastqs[i], read)
                sample.reads.append(read)

        all_samples.append(sample)
    return all_samples


def expand_metadata_per_sequence(
    df: pd.DataFrame,
    assembly_fasta: Path,
) -> pd.DataFrame:
    """
    Returns the sample's metadata row repeated once per FASTA record/header, with the record ID as the sequence name.
    seqsender matches the required `sequence_name` column to the FASTA record ID. Overwrites any existing
    `sequence_name`/`gb-sample_name` columns previously stored in the input metadata/terra table.
    """
    records = list(SeqIO.parse(assembly_fasta, "fasta"))
    if not records or any(len(record.seq) == 0 for record in records):
        raise ValueError("assembly FASTA is empty or has an empty sequence")
    seq_ids = [record.id for record in records]
    logger.info(f"Found {len(seq_ids)} FASTA record IDs in {assembly_fasta}")

    # one copy of the sample's row per sequence
    sample_metadata_df = pd.concat([df] * len(seq_ids), ignore_index=True)
    sample_metadata_df["sequence_name"] = seq_ids
    sample_metadata_df["gb-sample_name"] = seq_ids
    logger.info("Setting `sequence_name` and `gb-sample_name` metadata columns for each FASTA record ID")

    return sample_metadata_df


def check_existing_accessions(
    df: pd.DataFrame,
    biosample: bool,
    sra: bool,
    genbank: bool,
) -> None:
    """
    Raises if the sample already has an accession for a selected database, since submitting it again would make a duplicate.
    Any `genbank_accession` counts, so a flu sample with only some segments accessioned also fails.
    """
    columns = []
    flags = []
    for selected, accession_column, flag in (
        (biosample, "biosample_accession", "submit_to_biosample"),
        (sra, "sra_accession", "submit_to_sra"),
        (genbank, "genbank_accession", "submit_to_genbank"),
    ):
        if selected and accession_column in df and df[accession_column].iloc[0]:
            columns.append(accession_column)
            flags.append(flag)
    if columns:
        raise ValueError(f"already has {', '.join(columns)}; turn off {', '.join(flags)}")


def add_existing_accessions(
    df: pd.DataFrame,
    biosample: bool,
    sra: bool,
    definition_line_accessions: bool,
) -> None:
    """
    Adds the BioSample/SRA accessions from an earlier run to a GenBank submission's metadata, matching what seqsender
    adds itself when all three databases are submitted together. Mimics seqsender's `Link_Sample_Between_NCBI_Databases` and
    `Add_Definition_Line_Accessions` configuration options when submitting just to GenBank (with prior successful submission to BioSample and SRA).
    Necessary for attempting to resubmit a failed Genbank submission when everything else has successfully processed.
    """
    tags = []
    for selected, accession_column, modifier in (
        (biosample, "biosample_accession", "BioSample"),
        (sra, "sra_accession", "SRA"),
    ):
        accession = df[accession_column].iloc[0] if accession_column in df else ""
        if not selected and accession:
            df[f"src-{modifier}"] = accession
            tags.append(f"[{modifier}={accession}]")
            logger.info(f"Setting `src-{modifier}` metadata column to: '{accession}'")

    if definition_line_accessions and tags:
        column = "gb-fasta_definition_line_modifiers"
        existing = df[column].iloc[0] if column in df else ""
        # added after any modifiers already defined in the metadata table (mimics seqsender)
        modifiers = " ".join([*existing.split(), *tags])
        df[column] = modifiers
        logger.info(f"Setting `{column}` metadata column to: '{modifiers}'")


def build_sample_metadata(
    sample: Sample,
    metadata_df: pd.DataFrame,
    sample_column: str,
    biosample: bool,
    sra: bool,
    genbank: bool,
    link_databases: bool,
    definition_line_accessions: bool,
    bioproject: str,
    read_bucket: Optional[str] = None,
) -> pd.DataFrame:
    """
    Returns the sample's row from the input Terra table (`--metadata_file`) with the columns build_and_prep automatically fills in itself
    """
    # sample_column is the Terra sample identifier: f"entity:{args.metadata_file.stem}_id"
    sample_metadata_df = metadata_df[metadata_df[sample_column] == sample.sample_id].copy()
    logger.info(f"Filling in metadata columns for {sample.sample_id}")

    # adds required metadata columns for ALL submissions
    sample_metadata_df["bioproject"] = bioproject
    logger.info(f"Setting `bioproject` metadata column to: '{bioproject}'")

    # adds required metadata columns for SRA submissions
    if sample.reads:
        location = "cloud" if read_bucket else "local"
        sample_metadata_df["sra-file_location"] = location
        logger.info(f"Setting `sra-file_location` metadata column to: '{location}'")

        for n, read in enumerate(sample.reads, start=1):
            read_path = f"{read_bucket}/{read.name}" if read_bucket else read.name
            sample_metadata_df[f"sra-file_{n}"] = read_path
            logger.info(f"Setting `sra-file_{n}` metadata column to: '{read_path}'")

    # links a GenBank submission to the BioSample/SRA accessioned in an earlier run (if applicable)
    if genbank and link_databases:
        add_existing_accessions(
            df=sample_metadata_df,
            biosample=biosample,
            sra=sra,
            definition_line_accessions=definition_line_accessions,
        )

    return sample_metadata_df


def write_sample_metadata(
    sample: Sample,
    organism: str,
    biosample: bool,
    sra: bool,
    genbank: bool,
) -> None:
    """
    Creates a 'builds/<sample>/<sample>.metadata.csv' metadata file from `sample.metadata`, one row per FASTA record ID (genbank only).
    Raises ValueError if the sample can't be submitted.
    """
    sample_metadata_df = sample.metadata

    # '<description>' text left over from metadata_template.tsv would otherwise be submitted to NCBI as the value
    placeholder_columns = [column for column in sample_metadata_df if sample_metadata_df[column].str.fullmatch(r"<.*>").any()]
    if placeholder_columns:
        raise ValueError(f"fill in or clear the metadata template's '<description>' text in: {', '.join(placeholder_columns)}")

    check_existing_accessions(
        df=sample_metadata_df,
        biosample=biosample,
        sra=sra,
        genbank=genbank,
    )

    # GenBank submissions for flu requires explicit organism name. Fail here rather than during submission.
    if genbank and organism == "FLU":
        sample_organism = sample_metadata_df["organism"].iloc[0] if "organism" in sample_metadata_df else ""
        if sample_organism.strip().lower() not in ("influenza a virus", "influenza b virus", "influenza c virus"):
            raise ValueError(f"Organism '{sample_organism}' isn't accepted for GenBank flu submissions; Organism must be one of: 'Influenza A virus', 'Influenza B virus' or 'Influenza C virus'")

    logger.info(f"Building metadata file: {sample.metadata_csv}")

    # adds required metadata columns AND rows for Genbank submissions
    if genbank:
        sample_metadata_df = expand_metadata_per_sequence(
            df=sample_metadata_df,
            assembly_fasta=sample.fasta
        )

    sample_metadata_df.to_csv(sample.metadata_csv, index=False)


def write_metadata_template(
    metadata_df: pd.DataFrame,
    sample_column: str,
    organism: str,
    biosample_package: Optional[str],
    biosample: bool,
    sra: bool,
    genbank: bool,
    template_file: Path,
) -> None:
    """
    Writes '<table_name>_template.tsv'. Contains the selected samples' required metadata columns to fill in and upload to the Terra table.
    Mimics the behavior of the SeqSender submission wizard; the columns come from the schemas outlined in ('shiny/templates/'). The `metadata_df`
    holds all existing/built sample metadata (including Terra values); empty cells hold '<column description>'.
    """
    templates_dir = Path(shutil.which("seqsender.py")).parent / "shiny" / "templates"
    template_files = ["config.seqsender.seqsender.schema_template.csv"]
    if biosample:
        template_files.append(f"config.biosample.{biosample_package}_template.csv")
    if sra:
        template_files.append("config.sra.sra.schema_template.csv")
    if genbank:
        # genbank.schema_template.csv is left out: it only requires sequence_name and gb-sample_name, which come from the FASTA record IDs
        template_files.append("config.genbank.genbank.flu.src.schema_template.csv" if organism == "FLU" else "config.genbank.genbank.src.schema_template.csv")

    missing = [name for name in template_files if not (templates_dir / name).exists()]
    if missing:
        raise SystemExit(f"ERROR: seqsender has no metadata template {', '.join(missing)}; check BioSample_Package in the config file")

    columns_df = pd.concat([pd.read_csv(templates_dir / name, dtype=str, keep_default_na=False) for name in template_files])
    columns_df = columns_df[columns_df["required_column"] != "Optional"]
    # build_sample_metadata fills these in from the read inputs
    columns_df = columns_df[~columns_df["column_name"].isin(["sra-file_location", "sra-file_1"])]
    columns_df = columns_df.drop_duplicates("column_name")

    # descriptions as one line of plain text (no HTML tags or quotes) so the TSV needs no quoting
    placeholders = {
        row.column_name: "<" + " ".join(re.sub(r"<[^>]+>", "", row.description).replace('"', "").split()) + ">"
        for row in columns_df.itertuples()
    }

    template_df = metadata_df.reindex(columns=[sample_column, *placeholders]).fillna("")
    for column, placeholder in placeholders.items():
        template_df[column] = template_df[column].mask(template_df[column].str.strip() == "", placeholder)

    template_df.to_csv(template_file, sep="\t", index=False)
    logger.info(f"Wrote {len(placeholders)} required metadata columns for {len(template_df)} samples to {template_file}")


def summarize_prep_errors(prep_log: Path) -> str:
    """Returns the distinct errors in a seqsender prep log, joined with '; '"""
    lines = [line.strip() for line in prep_log.read_text().splitlines() if line.strip()]
    # flu samples have one row per segment, so the same error repeats with only the row index changed
    errors = [
        re.sub(r" at index '\d+'", "", line.removeprefix("Error: "))
        for line in lines
        if line.startswith("Error") and not line.startswith("Error: file ")
    ]
    if errors:
        return "; ".join(dict.fromkeys(errors))
    return lines[-1] if lines else "seqsender prep failed without output"


def run_prep(
    sample: Sample,
    db_flags: list[str],
    organism: str,
    config_file: Path,
) -> None:
    """Runs seqsender prep on the sample's folder and records whether it passed"""
    command = [
        "micromamba", "run", "-n", "seqsender", "seqsender.py", "prep", *db_flags,
        "--organism", organism,
        "--submission_dir", str(sample.sample_dir.parent),
        "--submission_name", sample.sample_id,
        "--config_file", str(config_file),
        "--metadata_file", str(sample.metadata_csv),
    ]
    if sample.fasta:
        command += ["--fasta_file", str(sample.fasta)]

    with open(sample.prep_log, "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)

    sample.passed = result.returncode == 0
    if not sample.passed:
        sample.reason = summarize_prep_errors(sample.prep_log)


def write_validation_summary(
    samples: list[Sample],
    output_dir: Path,
) -> None:
    """Writes every sample's PASS/FAIL status and reason"""
    summary_df = pd.DataFrame({
        "sample_id": [sample.sample_id for sample in samples],
        "status": ["PASS" if sample.passed else "FAIL" for sample in samples],
        "reason": [sample.reason for sample in samples],
    })
    summary_df.to_csv(output_dir / "validation_summary.tsv", sep="\t", index=False)

    num_passed = sum(sample.passed for sample in samples)
    (output_dir / "NUM_PASSED").write_text(f"{num_passed}")
    (output_dir / "NUM_FAILED").write_text(f"{len(samples) - num_passed}")
    logger.info(f"{num_passed} of {len(samples)} samples passed seqsender prep")


def upload_reads(
    valid_samples: list[Sample],
    read_bucket: str,
) -> None:
    """Copies the passing samples' reads to the transfer bucket NCBI fetches them from"""
    reads = [read for sample in valid_samples for read in sample.reads]
    logger.info(f"Copying {len(reads)} read files to {read_bucket}")

    bucket_name, _, prefix = read_bucket.removeprefix("gs://").partition("/")
    bucket = storage.Client().bucket(bucket_name)
    # FASTQs already in the `sra_transfer_bucket` with the same name are kept rather than overwritten;
    # any other upload error is raised
    transfer_manager.upload_many(
        [(str(read), bucket.blob(f"{prefix}/{read.name}".lstrip("/"))) for read in reads],
        skip_if_exists=True,
        raise_exception=True,
    )
    logger.info(f"Copied {len(reads)} read files to {read_bucket}")


def combine_passed(
    valid_samples: list[Sample],
    submit_dir: Path,
    config_file: Path,
    include_reads: bool = True,
) -> None:
    """
    Combines the passing samples' metadata, sequences and reads, plus the config, into 'submit/'.
    The reads are left out when NCBI fetches them from a transfer bucket (`include_reads=False`).
    """

    logger.info("Creating final submission folder for samples that passed seqsender prep validation")

    raw_reads_dir = submit_dir / "raw_reads"
    if include_reads:
        raw_reads_dir.mkdir(parents=True, exist_ok=True)

    combined_metadata_file = submit_dir / "metadata.csv"
    combined_fasta_file = submit_dir / "sequence.fasta"

    all_metadata: list[pd.DataFrame] = []
    all_fastas: list[SeqRecord] = []
    for sample in valid_samples:
        if include_reads:
            for read in sample.reads:
                logger.info(f"Moving {read} to {raw_reads_dir / read.name}")
                shutil.move(read, raw_reads_dir / read.name)
        all_metadata.append(pd.read_csv(sample.metadata_csv, dtype=str, keep_default_na=False))
        if sample.fasta:
            all_fastas.extend(SeqIO.parse(sample.fasta, "fasta"))

    combined_df = pd.concat(all_metadata, ignore_index=True)

    # moving these headers to the front for easier debugging
    front_cols = [c for c in ["bs-sample_name", "sra-sample_name", "gb-sample_name", "sequence_name"] if c in combined_df.columns]
    for i, col in enumerate(front_cols, start=1):
        combined_df.insert(i, col, combined_df.pop(col))

    logger.info(f"Combined {len(all_metadata)} sample metadata files into: {combined_metadata_file}")
    combined_df.to_csv(combined_metadata_file, index=False)

    if all_fastas:
        SeqIO.write(all_fastas, combined_fasta_file, "fasta")
        logger.info(f"Combined {len(all_fastas)} FASTA record IDs into: {combined_fasta_file}")

    logger.info(f"Copying {config_file} to {submit_dir / 'config.yaml'}")
    shutil.copy(config_file, submit_dir / "config.yaml")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output_dir / "seqsender_build_and_prep.log")
    db_flags = check_inputs(args)
    logger.info(f"Preparing {len(args.sample_ids)} samples for {', '.join(flag.lstrip('-') for flag in db_flags)}")

    builds_dir = args.output_dir / "builds"
    builds_dir.mkdir(parents=True, exist_ok=True)
    submit_dir = args.output_dir / "submit"
    submit_dir.mkdir(parents=True, exist_ok=True)

    # check the input config file for options that will modify our metadata/assembly fasta submission files
    ncbi_config = yaml.safe_load(args.config_file.read_text())["Submission"]["NCBI"]
    link_databases = ncbi_config["Link_Sample_Between_NCBI_Databases"]
    definition_line_accessions = ncbi_config.get("Add_Definition_Line_Accessions")

    # 1. stage only the files the selected databases use
    samples = stage_samples(
        builds_dir=builds_dir,
        sample_ids=args.sample_ids,
        assembly_fastas=args.assembly_fastas if args.genbank else [],
        read1_fastqs=args.read1_fastqs if args.sra else [],
        read2_fastqs=args.read2_fastqs if args.sra else [],
    )
    num_fastas = sum(sample.fasta is not None for sample in samples)
    num_reads = sum(len(sample.reads) for sample in samples)
    logger.info(f"Staged {len(samples)} samples ({num_fastas} FASTAs, {num_reads} reads) in {builds_dir}")

    # 2. the Terra table is the base of each sample's metadata.csv; entity:{table_name}_id holds the sample IDs
    metadata_df = pd.read_csv(args.metadata_file, sep="\t", dtype=str, keep_default_na=False)
    sample_column = f"entity:{args.metadata_file.stem}_id"

    # Excel re-saves 'YYYY-MM-DD' dates as 'M/D/YYYY', which seqsender rejects; any other value is left for seqsender prep to validate
    if "collection_date" in metadata_df:
        dates = pd.to_datetime(metadata_df["collection_date"], format="%m/%d/%Y", errors="coerce")
        metadata_df.loc[dates.notna(), "collection_date"] = dates.dt.strftime("%Y-%m-%d")
        logger.info(f"Converted {dates.notna().sum()} M/D/YYYY `collection_date` values to YYYY-MM-DD")

    # 3. fill in the columns build_and_prep sets itself; the samples are checked in step 5
    for sample in samples:
        sample.metadata = build_sample_metadata(
            sample=sample,
            metadata_df=metadata_df,
            sample_column=sample_column,
            biosample=args.biosample,
            sra=args.sra,
            genbank=args.genbank,
            link_databases=link_databases,
            definition_line_accessions=definition_line_accessions,
            bioproject=args.bioproject,
            read_bucket=args.sra_transfer_bucket,
        )

    # 4. the template of required columns is written even if every sample fails validation below
    write_metadata_template(
        metadata_df=pd.concat([sample.metadata for sample in samples]),
        sample_column=sample_column,
        organism=args.organism,
        biosample_package=ncbi_config.get("BioSample_Package"),
        biosample=args.biosample,
        sra=args.sra,
        genbank=args.genbank,
        template_file=args.output_dir / f"{args.metadata_file.stem}_template.tsv",
    )

    # 5. validate each sample on its own; a sample that fails is left out of the submission
    for sample in samples:
        try:
            write_sample_metadata(
                sample=sample,
                organism=args.organism,
                biosample=args.biosample,
                sra=args.sra,
                genbank=args.genbank,
            )
        except ValueError as error:
            sample.reason = str(error)
        else:
            run_prep(
                sample=sample,
                db_flags=db_flags,
                organism=args.organism,
                config_file=args.config_file,
            )

        if sample.passed:
            logger.info(f"{sample.sample_id}: PASS")
        else:
            logger.warning(f"{sample.sample_id}: FAIL: {sample.reason}")

    write_validation_summary(samples, args.output_dir)

    # 6. copy valid samples to the `--sra_transfer_bucket` if appropriate. exit if no valid samples
    valid_samples = [sample for sample in samples if sample.passed]
    if not valid_samples:
        logger.warning("No samples passed seqsender prep; no submission tarball written")
        return

    if args.sra and args.sra_transfer_bucket:
        if args.validation_only:
            logger.info(f"--validation_only: the reads were not copied to {args.sra_transfer_bucket}")
        else:
            upload_reads(
                valid_samples=valid_samples,
                read_bucket=args.sra_transfer_bucket,
            )

    # 7. combine the passing samples into submit/ and tar it for the submit task
    combine_passed(
        valid_samples=valid_samples,
        submit_dir=submit_dir,
        config_file=args.config_file,
        include_reads=not args.sra_transfer_bucket,
    )

    # reads are already gzipped, so compress lightly
    tarball = args.output_dir / "seqsender_post_validation_checkpoint.tar.gz"
    with tarfile.open(tarball, "w:gz", compresslevel=1) as tar:
        tar.add(submit_dir, arcname="submit")
    logger.info(f"Wrote {tarball}")


if __name__ == "__main__":
    main()
