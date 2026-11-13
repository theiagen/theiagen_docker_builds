import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .models import BIOSAMPLE_DB, SRA_DB, BuildOptions, Sample, SampleError

logger = logging.getLogger(__name__)


def get_value(
    df: pd.DataFrame,
    column: str,
) -> str:
    """Returns the sample's value in `column`, or '' if the column isn't in its metadata"""
    return df[column].iloc[0] if column in df else ""


def expand_metadata_per_sequence(
    df: pd.DataFrame,
    assembly_fasta: Path,
) -> pd.DataFrame:
    """
    Returns the sample's metadata row repeated once per FASTA record/header, with the record ID as the sequence name.
    seqsender matches the required `sequence_name` column to the FASTA record ID. Overwrites any existing
    `sequence_name`/`gb-sample_name` columns previously stored in the input metadata/terra table.
    """
    try:
        records = list(SeqIO.parse(assembly_fasta, "fasta"))
    except ValueError as error:
        raise SampleError(f"assembly FASTA can't be read: {error}") from error
    if not records or any(len(record.seq) == 0 for record in records):
        raise SampleError("assembly FASTA is empty or has an empty sequence")
    seq_ids = [record.id for record in records]
    logger.info(f"Found {len(seq_ids)} FASTA record IDs in {assembly_fasta}")

    # one copy of the sample's row per sequence
    sample_metadata_df = pd.concat([df] * len(seq_ids), ignore_index=True)
    sample_metadata_df["sequence_name"] = seq_ids
    sample_metadata_df["gb-sample_name"] = seq_ids
    logger.info("Setting `sequence_name` and `gb-sample_name` metadata columns for each FASTA record ID")

    return sample_metadata_df


def add_existing_accessions(
    df: pd.DataFrame,
    options: BuildOptions,
) -> None:
    """
    Adds the BioSample/SRA accessions from an earlier run to a GenBank submission's metadata, matching what seqsender
    adds itself when all three databases are submitted together. Mimics seqsender's `Link_Sample_Between_NCBI_Databases` and
    `Add_Definition_Line_Accessions` configuration options when submitting just to GenBank (with prior successful submission to BioSample and SRA).
    Necessary for attempting to resubmit a failed Genbank submission when everything else has successfully processed.
    """
    tags = []

    biosample_accession = get_value(df, BIOSAMPLE_DB.accession_column)
    # if `biosample_accession` column exists, but `submit_to_biosample = false`; we can still link existing accessions to GenBank records
    if biosample_accession and not options.biosample:
        df["src-BioSample"] = biosample_accession
        tags.append(f"[BioSample={biosample_accession}]")
        logger.info(f"Setting `src-BioSample` metadata column to: '{biosample_accession}'")

    # if `sra_accession` column exists, but `submit_to_sra = false`; we can still link existing accessions to GenBank records
    sra_accession = get_value(df, SRA_DB.accession_column)
    if sra_accession and not options.sra:
        df["src-SRA"] = sra_accession
        tags.append(f"[SRA={sra_accession}]")
        logger.info(f"Setting `src-SRA` metadata column to: '{sra_accession}'")

    if options.definition_line_accessions and tags:
        # if existing `gb-fasta_definition_line_modifiers` column; use that
        def_line_accession = get_value(df, "gb-fasta_definition_line_modifiers")
        # added after any modifiers already defined in the metadata table (mimics seqsender)
        modifiers = " ".join([*def_line_accession.split(), *tags])
        df["gb-fasta_definition_line_modifiers"] = modifiers
        logger.info(f"Setting `gb-fasta_definition_line_modifiers` metadata column to: '{modifiers}'")


def build_sample_metadata(
    sample: Sample,
    metadata_df: pd.DataFrame,
    sample_column: str,
    options: BuildOptions,
) -> pd.DataFrame:
    """
    Returns the sample's row from the input Terra table (`--metadata_file`) with the columns build_and_prep automatically fills in itself
    """
    # sample_column is the Terra sample identifier: f"entity:{args.metadata_file.stem}_id"
    sample_metadata_df = metadata_df[metadata_df[sample_column] == sample.sample_id].copy()
    logger.info(f"Filling in metadata columns for {sample.sample_id}")

    # Excel re-saves 'YYYY-MM-DD' dates as 'M/D/YYYY', which seqsender rejects; any other value is left for seqsender prep to validate
    if "collection_date" in sample_metadata_df:
        dates = pd.to_datetime(sample_metadata_df["collection_date"], format="%m/%d/%Y", errors="coerce")
        if dates.notna().any():
            converted = dates.dt.strftime("%Y-%m-%d")
            sample_metadata_df.loc[dates.notna(), "collection_date"] = converted
            logger.info(f"Setting `collection_date` metadata column to: '{converted.iloc[0]}' (converted from M/D/YYYY)")

    # adds required metadata columns for ALL submissions
    sample_metadata_df["bioproject"] = options.bioproject
    logger.info(f"Setting `bioproject` metadata column to: '{options.bioproject}'")

    # adds required metadata columns for SRA submissions
    if sample.reads:
        location = "cloud" if options.sra_transfer_bucket else "local"
        sample_metadata_df["sra-file_location"] = location
        logger.info(f"Setting `sra-file_location` metadata column to: '{location}'")

        for n, read in enumerate(sample.reads, start=1):
            read_path = f"{options.sra_transfer_bucket}/{read.name}" if options.sra_transfer_bucket else read.name
            sample_metadata_df[f"sra-file_{n}"] = read_path
            logger.info(f"Setting `sra-file_{n}` metadata column to: '{read_path}'")

    # links a GenBank submission to the BioSample/SRA accessioned in an earlier run (if applicable)
    if options.genbank and options.link_databases:
        add_existing_accessions(
            df=sample_metadata_df,
            options=options,
        )

    return sample_metadata_df


def write_sample_metadata(
    sample: Sample,
    options: BuildOptions,
) -> None:
    """
    Creates a 'builds/<sample>/<sample>.metadata.csv' metadata file from `sample.metadata`, one row per FASTA record ID (genbank only).
    Raises SampleError if the sample's assembly FASTA is empty or can't be read.
    """
    sample_metadata_df = sample.metadata

    logger.info(f"Building metadata file: {sample.metadata_csv}")

    # adds required metadata columns AND rows for Genbank submissions
    if options.genbank:
        sample_metadata_df = expand_metadata_per_sequence(
            df=sample_metadata_df,
            assembly_fasta=sample.fasta
        )

    sample_metadata_df.to_csv(sample.metadata_csv, index=False)
