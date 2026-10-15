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
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
    parser.add_argument("--sample_ids", nargs="+", required=True, help="Sample IDs to be submitted")
    parser.add_argument("--metadata_file", type=Path, required=True, help="Terra table TSV, named <table_name>.tsv")
    parser.add_argument("--organism", type=str, required=True, help="seqsender organism, e.g. FLU")
    parser.add_argument("--config_file", type=Path, required=True, help="Path to the seqsender configuration file")

    # Terra lists these in the same order as the sample IDs
    parser.add_argument("--assembly_fastas", nargs="*", type=Path, default=[], help="Required for GenBank; one assembly FASTA per sample")
    parser.add_argument("--read1_fastqs", nargs="*", type=Path, default=[], help="Required for SRA; one 'read1' FASTQ per sample")
    parser.add_argument("--read2_fastqs", nargs="*", type=Path, default=[], help="Optional for SRA; leave out for single-end reads; one 'read2' FASTQ per sample")

    parser.add_argument("--biosample", action="store_true")
    parser.add_argument("--sra", action="store_true")
    parser.add_argument("--genbank", action="store_true")

    parser.add_argument("--output_dir", type=Path, default=Path.cwd(), help="Directory to write the submission files to; populates 'builds' and 'submit' directories")
    return parser.parse_args()


@dataclass
class Sample:
    """A sample's staged files and its seqsender prep result"""
    sample_id: str
    sample_dir: Path  # builds/<sample_id>
    fasta: Optional[Path] = None  # builds/<sample_id>/<sample_id>.assembly.fasta
    reads: list[Path] = field(default_factory=list)  # in builds/raw_reads, read 1 first
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


def write_sample_metadata(
    sample: Sample,
    metadata_df: pd.DataFrame,
    sample_column: str,
    genbank: bool
) -> None:
    """
    Creates a 'builds/<sample>/<sample>.metadata.csv' metadata file for each input `sample_id`.
    Each metadata.csv file is based on the original input Terra table (`--metadata_file`), one row per FASTA record ID (genbank only).
    Fills in required columns for submission depending on the selected databases.
    """
    # sample_column is the Terra sample identifier: f"entity:{args.metadata_file.stem}_id"
    sample_metadata_df = metadata_df[metadata_df[sample_column] == sample.sample_id].copy()

    if len(sample_metadata_df) != 1:
        raise ValueError(f"found {len(sample_metadata_df)} rows for this sample in the table")

    logger.info(f"Building metadata file: {sample.metadata_csv}")

    # adds required metadata columns for SRA submissions
    if sample.reads:
        sample_metadata_df["sra-file_location"] = "local"
        logger.info("Setting `sra-file_location` metadata column to: 'local'")

        for n, read in enumerate(sample.reads, start=1):
            sample_metadata_df[f"sra-file_{n}"] = read.name
            logger.info(f"Setting `sra-file_{n}` metadata column to: '{read.name}'")

    # adds required metadata columns AND rows for Genbank submissions
    if genbank:
        sample_metadata_df = expand_metadata_per_sequence(
            df=sample_metadata_df,
            assembly_fasta=sample.fasta
        )

    sample_metadata_df.to_csv(sample.metadata_csv, index=False)


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
    config_file: Path
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
    output_dir: Path
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


def combine_passed(
    passed: list[Sample],
    submit_dir: Path,
    config_file: Path
) -> None:
    """Combines the passing samples' metadata, sequences and reads, plus the config, into 'submit/'"""

    logger.info("Creating final submission folder for samples that passed seqsender prep validation")

    raw_reads_dir = submit_dir / "raw_reads"
    raw_reads_dir.mkdir(parents=True, exist_ok=True)

    combined_metadata_file = submit_dir / "metadata.csv"
    combined_fasta_file = submit_dir / "sequence.fasta"

    all_metadata: list[pd.DataFrame] = []
    all_fastas: list[SeqRecord] = []
    for sample in passed:
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
    submit_dir = args.output_dir / "submit"

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

    # 3. validate each sample on its own; a sample that fails is left out of the submission
    for sample in samples:
        try:
            write_sample_metadata(
                sample=sample,
                metadata_df=metadata_df,
                sample_column=sample_column,
                genbank=args.genbank
            )
        except ValueError as error:
            sample.reason = str(error)
        else:
            run_prep(
                sample=sample,
                db_flags=db_flags,
                organism=args.organism,
                config_file=args.config_file
            )

        if sample.passed:
            logger.info(f"{sample.sample_id}: PASS")
        else:
            logger.warning(f"{sample.sample_id}: FAIL: {sample.reason}")

    write_validation_summary(samples, args.output_dir)

    # 4. combine the passing samples into submit/ and tar it for the submit task
    passed = [sample for sample in samples if sample.passed]
    if not passed:
        logger.warning("No samples passed seqsender prep; no submission tarball written")
        return
    combine_passed(passed, submit_dir, args.config_file)

    # reads are already gzipped, so compress lightly
    tarball = args.output_dir / "seqsender_build_and_prep.tar.gz"
    with tarfile.open(tarball, "w:gz", compresslevel=1) as tar:
        tar.add(submit_dir, arcname="submit")
    logger.info(f"Wrote {tarball}")


if __name__ == "__main__":
    main()
