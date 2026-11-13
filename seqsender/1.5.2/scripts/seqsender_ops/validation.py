import logging
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from .metadata import get_value
from .models import BuildOptions, Sample, SampleError

logger = logging.getLogger(__name__)


def check_existing_accessions(
    df: pd.DataFrame,
    options: BuildOptions,
) -> None:
    """
    Raises if the sample already has an accession for a selected database, since submitting it again would make a duplicate.
    Any value(s) under the `genbank_accession` column will count. This allows flu samples with partial accessions assigned to also fail.
    """
    # <db>_accession column exists and is non-empty
    valid_accessions = [db for db in options.target_dbs if get_value(df, db.accession_column)]
    if valid_accessions:
        raise SampleError(
            f"already has {', '.join(db.accession_column for db in valid_accessions)}; "
            f"turn off {', '.join(f'submit_to_{db.name}' for db in valid_accessions)}"
        )


def validate_sample_metadata(
    sample_metadata_df: pd.DataFrame,
    options: BuildOptions,
) -> None:
    """
    Checks the sample's metadata for problems `seqsender prep` doesn't catch, so they fail here rather than during submission.
    Raises SampleError if the sample can't be submitted.
    """
    # '<description>' text left over from metadata_template.tsv would otherwise be submitted to NCBI as the value
    placeholder_columns = [column for column in sample_metadata_df if sample_metadata_df[column].str.fullmatch(r"<.*>").any()]
    if placeholder_columns:
        raise SampleError(f"fill in or clear the metadata template's '<description>' text in: {', '.join(placeholder_columns)}")

    check_existing_accessions(
        df=sample_metadata_df,
        options=options,
    )

    # GenBank submissions for flu requires explicit organism name. Fail here rather than during submission.
    if options.genbank and options.seqsender_organism_tag == "FLU":
        sample_organism = get_value(sample_metadata_df, "organism")
        if sample_organism.strip().lower() not in ("influenza a virus", "influenza b virus", "influenza c virus"):
            raise SampleError(f"Organism '{sample_organism}' isn't accepted for GenBank flu submissions; Organism must be one of: 'Influenza A virus', 'Influenza B virus' or 'Influenza C virus'")


def write_metadata_template(
    metadata_df: pd.DataFrame,
    sample_column: str,
    options: BuildOptions,
    template_file: Path,
) -> None:
    """
    Writes '<table_name>_template.tsv'. Contains the selected samples' required metadata columns to fill in and upload to the Terra table.
    Mimics the behavior of the SeqSender submission wizard; the columns come from the schemas outlined in ('shiny/templates/'). The `metadata_df`
    holds all existing/built sample metadata (including Terra values); empty cells hold '<column description>'.
    """
    templates_dir = Path(shutil.which("seqsender.py")).parent / "shiny" / "templates"
    template_files = ["config.seqsender.seqsender.schema_template.csv"]
    if options.biosample:
        template_files.append(f"config.biosample.{options.biosample_package}_template.csv")
    if options.sra:
        template_files.append("config.sra.sra.schema_template.csv")
    if options.genbank:
        # genbank.schema_template.csv is left out: it only requires sequence_name and gb-sample_name, which come from the FASTA record IDs
        template_files.append("config.genbank.genbank.flu.src.schema_template.csv" if options.seqsender_organism_tag == "FLU" else "config.genbank.genbank.src.schema_template.csv")

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
    options: BuildOptions,
) -> None:
    """Runs seqsender prep on the sample's folder; raises SampleError with seqsender's errors if it fails"""
    command = [
        "micromamba", "run", "-n", "seqsender", "seqsender.py", "prep", *options.db_flags,
        "--organism", options.seqsender_organism_tag,
        "--submission_dir", str(sample.sample_dir.parent),
        "--submission_name", sample.sample_id,
        "--config_file", str(options.config_file),
        "--metadata_file", str(sample.metadata_csv),
    ]
    if sample.fasta:
        command += ["--fasta_file", str(sample.fasta)]

    with open(sample.prep_log, "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        raise SampleError(summarize_prep_errors(sample.prep_log))


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
