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
