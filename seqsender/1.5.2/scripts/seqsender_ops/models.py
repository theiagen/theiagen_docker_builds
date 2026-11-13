from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd


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


@dataclass(frozen=True)
class Database:
    """
    Immutable identifiers associated with each submission database.
    Used across seqsender, metadata columns, WDL workflow and Terra.
    """
    name: str  # seqsender CLI flag, without the '--': "biosample"
    sample_name_column: str  # submission_status_report.csv: "bs-sample_name"
    accession_column: str  # submission_status_report.csv: "biosample_accession"
    terra_toggle: str  # WDL input that selects the database: "submit_to_biosample"
    src_modifier: Optional[str]  # GenBank source modifier linking a sample to an earlier accession: "src-BioSample"


BIOSAMPLE_DB = Database(
    name="biosample",
    sample_name_column="bs-sample_name",
    accession_column="biosample_accession",
    terra_toggle="submit_to_biosample",
    src_modifier="BioSample"
)

SRA_DB = Database(
    name="sra",
    sample_name_column="sra-sample_name",
    accession_column="sra_accession",
    terra_toggle="submit_to_sra",
    src_modifier="SRA"
)

GENBANK_DB = Database(
    name="genbank",
    sample_name_column="gb-sample_name",
    accession_column="genbank_accession",
    terra_toggle="submit_to_genbank",
    src_modifier=None
)

DATABASES = (BIOSAMPLE_DB, SRA_DB, GENBANK_DB)
