import logging
import re
import subprocess
from pathlib import Path

import pandas as pd

from .metadata import flu_genbank_required_columns, get_value
from .models import BuildOptions, Sample, SampleError

logger = logging.getLogger(__name__)


def check_placeholders(
    df: pd.DataFrame,
) -> None:
    """Raises if '<description>' text from the metadata template is left in a column, since it would be submitted to NCBI as the value"""
    placeholder_columns = [column for column in df if df[column].str.fullmatch(r"<.*>").any()]
    if placeholder_columns:
        raise SampleError(f"fill in or clear the metadata template's '<description>' text in: {', '.join(placeholder_columns)}")


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


def check_flu_genbank(
    df: pd.DataFrame,
) -> None:
    """
    Raises if a GenBank flu sample is missing what NCBI requires and `seqsender prep` doesn't check:
    an accepted organism name and the columns from `flu_genbank_required_columns`
    """
    sample_organism = get_value(df, "organism")
    if sample_organism.strip().lower() not in ("influenza a virus", "influenza b virus", "influenza c virus"):
        raise SampleError(f"Organism '{sample_organism}' isn't accepted for GenBank flu submissions; Organism must be one of: 'Influenza A virus', 'Influenza B virus' or 'Influenza C virus'")

    missing_flu_columns = [column for column in flu_genbank_required_columns(sample_organism) if not get_value(df, column).strip()]
    if missing_flu_columns:
        raise SampleError(f"NCBI requires {', '.join(missing_flu_columns)} for '{sample_organism}' submissions to GenBank")


def validate_sample_metadata(
    sample_metadata_df: pd.DataFrame,
    options: BuildOptions,
) -> None:
    """
    Checks the sample's metadata for problems `seqsender prep` doesn't catch, so they fail here rather than during submission.
    Raises SampleError if the sample can't be submitted.
    """
    check_placeholders(
        df=sample_metadata_df,
    )
    check_existing_accessions(
        df=sample_metadata_df,
        options=options,
    )
    if options.genbank and options.seqsender_organism_tag == "FLU":
        check_flu_genbank(
            df=sample_metadata_df,
        )


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
        "reason": [sample.error or "" for sample in samples],
    })
    summary_df.to_csv(output_dir / "validation_summary.tsv", sep="\t", index=False)

    num_passed = sum(sample.passed for sample in samples)
    (output_dir / "NUM_PASSED").write_text(f"{num_passed}")
    (output_dir / "NUM_FAILED").write_text(f"{len(samples) - num_passed}")
    logger.info(f"{num_passed} of {len(samples)} samples passed seqsender prep")
