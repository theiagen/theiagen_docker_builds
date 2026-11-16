import logging
import re
import subprocess
from pathlib import Path
from typing import Optional

from .models import BuildOptions, SampleError

logger = logging.getLogger(__name__)


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


def seqsender_prep(
    submission_dir: Path,
    submission_name: str,
    metadata_file: Path,
    fasta_file: Optional[Path],
    prep_log: Path,
    options: BuildOptions,
) -> None:
    """
    Runs `seqsender prep`, which writes its files to '<submission_dir>/<submission_name>/', and writes its output to `prep_log`.
    Raises SampleError with seqsender's errors if it fails.
    """
    command = [
        "micromamba", "run", "-n", "seqsender", "seqsender.py", "prep", *options.db_flags,
        "--organism", options.seqsender_organism_tag,
        "--submission_dir", str(submission_dir),
        "--submission_name", submission_name,
        "--config_file", str(options.config_file),
        "--metadata_file", str(metadata_file),
    ]
    if fasta_file:
        command += ["--fasta_file", str(fasta_file)]

    with open(prep_log, "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        raise SampleError(summarize_prep_errors(prep_log))


def seqsender_submission_status(
    submit_dir: Path,
    submission_name: str,
    key: str,
) -> None:
    """Runs one `seqsender submission_status` check; exits if seqsender reports any problem"""
    command = [
        "micromamba", "run", "-n", "seqsender", "seqsender.py", "submission_status",
        "--submission_dir", str(submit_dir),
        "--submission_name", submission_name,
        "--key", key,
    ]
    # the command isn't logged because it holds the key
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    for line in lines:
        logger.info(f"seqsender: {line}")

    # seqsender prints these and can still exit 0, e.g. when it can't read the accessions from a processed report.xml
    if "Unable to access config file with key" in result.stdout:
        raise SystemExit("ERROR: --key doesn't decrypt the NCBI password in config.yaml; use the key that encrypted it")
    if result.returncode != 0 or any(line.lstrip().startswith("Error") for line in lines):
        raise SystemExit("ERROR: seqsender submission_status failed; see the seqsender output in the log")