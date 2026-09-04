from dataclasses import dataclass, field, fields

from ngs_pipeline_lib.tools.tools import hash_sequence


@dataclass
class AlleleCall:
    seq: str
    seq_up: str
    seq_down: str
    identity: float
    evalue: float
    bitscore: int
    align_len: int
    has_complete_alignment: bool
    repeat_score: float
    mismatch_count: int
    gap_count: int
    contig: str
    pos_start: int
    pos_end: int
    is_strand_forward: bool
    non_atgc_count: int
    has_codon_start: bool
    has_codon_stop: bool
    is_start_stop_codon_required: bool
    has_internal_codon_stop: bool

    id: int = field(init=False)
    flanking_seq: str = field(init=False)
    flanking_id: int = field(init=False)

    def __post_init__(self):
        self.id = hash_sequence(self.seq)
        self.flanking_seq = self.seq_up + self.seq + self.seq_down
        self.flanking_id = hash_sequence(self.flanking_seq)

        for field_ in fields(self):
            setattr(self, field_.name, field_.type(getattr(self, field_.name)))

    @property
    def has_codon_issues(self) -> bool:
        return (
            self.is_start_stop_codon_required
            and (not self.has_codon_start or not self.has_codon_stop)
        ) or self.has_internal_codon_stop


@dataclass
class FlaggedAlleleCall:
    id: int
    seq: str
    flag: int
    flank_id: int
    flank_seq: str
    coverage_prob: int | None
    bias_prob: int | None
    double_calls: int | None
    more_calls: int | None
    mismatch_calls: int | None
