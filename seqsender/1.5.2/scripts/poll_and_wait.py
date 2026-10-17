#!/usr/bin/env python3

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

logger = logging.getLogger("poll_and_wait")

# DELETED: NCBI removed the submission, e.g. it stayed in error too long
ERROR_STATUSES = {"ERROR", "FAILED", "DELETED"}
# submission_status_report.csv sample name and accession columns for each database
ACCESSION_COLUMNS = {
    "BIOSAMPLE": ("bs-sample_name", "biosample_accession"),
    "SRA": ("sra-sample_name", "sra_accession"),
    "GENBANK": ("gb-sample_name", "genbank_accession"),
}


def setup_logging(log_file: Path) -> None:
    """Logs to stdout and to log_file"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode="a")],
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Polls NCBI for a seqsender submission until it finishes or times out")
    parser.add_argument("--submit_dir", type=Path, required=True, help="The extracted 'submit' directory holding submission_log.csv")
    parser.add_argument("--key", type=str, required=True, help="Key that decrypts the NCBI password in the submission's config file")
    parser.add_argument("--max_wait_hours", type=float, default=24, help="Stops with PENDING after this many hours")
    parser.add_argument("--poll_interval_minutes", type=float, default=10, help="Minutes between checks")
    parser.add_argument("--output_dir", type=Path, default=Path.cwd(), help="Directory to write the log, FINAL_STATUS and reports to")
    return parser.parse_args()


def load_submission_log(submit_dir: Path) -> pd.DataFrame:
    """Reads submission_log.csv as strings, so it's written back unchanged"""
    return pd.read_csv(submit_dir / "submission_log.csv", dtype=str)


def relocate_submission_log(
    submit_dir: Path,
    log_df: pd.DataFrame,
) -> None:
    """Points submission_log.csv's absolute paths at submit_dir when the state was extracted somewhere else"""
    # seqsender records absolute paths when it submits; Submission_Directory is <submit_dir>/<submission_name>/submission_files/<DATABASE>
    logged_dir = str(Path(log_df["Submission_Directory"].iloc[0]).parents[2])
    if logged_dir != str(submit_dir):
        for column in ("Submission_Directory", "Config_File"):
            log_df[column] = log_df[column].str.replace(logged_dir, str(submit_dir), n=1, regex=False)
        log_df.to_csv(submit_dir / "submission_log.csv", index=False)
        logger.info(f"Moved the submission_log.csv paths from {logged_dir} to {submit_dir}")

    for path in (log_df["Submission_Directory"].iloc[0], log_df["Config_File"].iloc[0]):
        if not Path(path).exists():
            raise SystemExit(f"ERROR: {path} from submission_log.csv doesn't exist")


def get_db_status(log_df: pd.DataFrame) -> dict[str, str]:
    """Returns each database's latest status, e.g. {"BIOSAMPLE": "PROCESSED", "SRA": "PROCESSING", "GENBANK": "WAITING"}"""
    # GENBANK-FTP and GENBANK-TBL2ASN are both reported as GENBANK
    return {database.split("-")[0]: status for database, status in zip(log_df["Database"], log_df["Submission_Status"])}


def get_overall_status(statuses: dict[str, str]) -> str:
    """Returns ERROR if any database has an error, COMPLETE if all are processed, otherwise PENDING"""
    if any(status in ERROR_STATUSES for status in statuses.values()):
        return "ERROR"
    if all(status == "PROCESSED" for status in statuses.values()):
        return "COMPLETE"
    return "PENDING"


def run_submission_status(
    submit_dir: Path,
    submission_name: str,
    key: str,
) -> None:
    """Runs one `seqsender.py submission_status` check; exits if seqsender reports any problem"""
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


def missing_accessions(status_report: Path) -> dict[str, list[str]]:
    """Returns the sample names in submission_status_report.csv that have no accession, for each database"""
    report_df = pd.read_csv(status_report, dtype=str, keep_default_na=False)
    missing = {}
    for database, (name_column, accession_column) in ACCESSION_COLUMNS.items():
        if name_column in report_df:
            names = report_df.loc[report_df[accession_column].str.strip() == "", name_column].unique().tolist()
            if names:
                missing[database] = names
    return missing


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output_dir / "seqsender_poll_and_wait.log")

    submit_dir = args.submit_dir.absolute()
    if not (submit_dir / "submission_log.csv").is_file():
        raise SystemExit(f"ERROR: {submit_dir}/submission_log.csv not found; this state hasn't been submitted")

    log_df = load_submission_log(submit_dir)
    submission_name = log_df["Submission_Name"].iloc[0]
    (args.output_dir / "SUBMISSION_NAME").write_text(submission_name)

    relocate_submission_log(
        submit_dir=submit_dir,
        log_df=log_df
    )

    logger.info(f"Polling {submission_name} ({', '.join(log_df['Database'])}) every {args.poll_interval_minutes:g} minutes for up to {args.max_wait_hours:g} hours")

    # the first check runs immediately; seqsender also uploads GenBank once BioSample and SRA are processed
    deadline = time.time() + args.max_wait_hours * 3600
    interval = args.poll_interval_minutes * 60
    while True:
        run_submission_status(
            submit_dir=submit_dir,
            submission_name=submission_name,
            key=args.key,
        )

        # seqsender updates submission_log.csv on each check
        log_df = load_submission_log(submit_dir)
        db_statuses = get_db_status(log_df)
        main_status = get_overall_status(db_statuses)
        logger.info(f"{' '.join(f'{db}={status}' for db, status in db_statuses.items())} -> {main_status}")

        if main_status != "PENDING" or time.time() + interval > deadline:
            break
        time.sleep(interval)

    # NCBI can report processed-ok without accessioning every sample, e.g. GenBank_Auto_Remove_Failed_Samples drops failed sequences
    if main_status == "COMPLETE":
        missing = missing_accessions(submit_dir / submission_name / "submission_files" / "submission_status_report.csv")
        for database, names in missing.items():
            logger.error(f"{database}: {len(names)} processed without an accession: {', '.join(names)}")
        if missing:
            main_status = "ERROR"

    (args.output_dir / "DATABASE_STATUSES").write_text(' '.join(f'{db}={status}' for db, status in db_statuses.items()))
    (args.output_dir / "FINAL_STATUS").write_text(main_status)
    logger.info(f"Final status: {main_status}")


if __name__ == "__main__":
    main()