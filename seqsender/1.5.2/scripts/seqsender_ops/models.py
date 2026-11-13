from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml


@dataclass(frozen=True)
class Database:
    """
    Immutable identifiers associated with each submission database.
    Used across seqsender, metadata columns, WDL workflow and Terra.
    """
    name: str  # seqsender CLI flag, without the '--': "biosample"
    sample_name_column: str  # submission_status_report.csv: "bs-sample_name"
    accession_column: str  # submission_status_report.csv: "biosample_accession"

BIOSAMPLE_DB = Database(
    name="biosample",
    sample_name_column="bs-sample_name",
    accession_column="biosample_accession",
)
SRA_DB = Database(
    name="sra",
    sample_name_column="sra-sample_name",
    accession_column="sra_accession",
)
GENBANK_DB = Database(
    name="genbank",
    sample_name_column="gb-sample_name",
    accession_column="genbank_accession",
)
DATABASES = (BIOSAMPLE_DB, SRA_DB, GENBANK_DB)


@dataclass(frozen=True)
class BuildOptions:
    """
    General run-wide settings for the `build_and_prep` script and seqsender workflow. Stores argparse settings from the command line
    and the input seqsender config file. Built in the main function and passed to each step instead of each individual flag.
    """
    config_file: Path
    bioproject: str
    biosample: bool
    sra: bool
    genbank: bool
    validation_only: bool
    seqsender_organism_tag: str  # passed to `seqsender prep --organism`: "FLU", "COV" or "OTHER"
    sra_transfer_bucket: Optional[str]  # 'gs://' path without a trailing '/'; None uploads the reads over FTP
    link_databases: bool  # config file `Link_Sample_Between_NCBI_Databases`
    definition_line_accessions: bool  # config file `Add_Definition_Line_Accessions`
    biosample_package: Optional[str]  # config file `BioSample_Package`

    @classmethod
    def from_args(cls, args) -> "BuildOptions":
        seqsender_config = yaml.safe_load(args.config_file.read_text())["Submission"]["NCBI"]
        return cls(
            config_file=args.config_file,
            bioproject=args.bioproject,
            biosample=args.biosample,
            sra=args.sra,
            genbank=args.genbank,
            validation_only=args.validation_only,
            seqsender_organism_tag=args.seqsender_organism_tag,
            sra_transfer_bucket=args.sra_transfer_bucket.rstrip("/") if args.sra_transfer_bucket else None,
            link_databases=seqsender_config["Link_Sample_Between_NCBI_Databases"],
            definition_line_accessions=bool(seqsender_config.get("Add_Definition_Line_Accessions")),
            biosample_package=seqsender_config.get("BioSample_Package"),
        )

    @property
    def target_dbs(self) -> tuple[Database, ...]:
        """The databases selected for this run, e.g. (BIOSAMPLE_DB, SRA_DB)"""
        return tuple(db for db in DATABASES if getattr(self, db.name))

    @property
    def db_flags(self) -> list[str]:
        """seqsender flags for the selected databases, e.g. ["--biosample", "--sra"]"""
        return [f"--{db.name}" for db in self.target_dbs]


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