import logging
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from google.cloud import storage
from google.cloud.storage import transfer_manager

from .models import Sample

logger = logging.getLogger(__name__)


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
