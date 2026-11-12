import logging
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO

from .models import Sample
from .validation import check_existing_accessions

logger = logging.getLogger(__name__)


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
    seqsender_organism_tag: str,
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
    if genbank and seqsender_organism_tag == "FLU":
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
