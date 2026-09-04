from enum import Enum
from typing import Self


class Flag(Enum):
    # Group
    CORE = (1, "Core locus", None)
    ACCESSORY = (2, "Accessory locus", None)
    # Called
    CALLED = (4, "Locus called", None)
    NOT_CALLED = (8, "Locus not called", None)
    NOT_FOUND = (16, "Locus not found", None)
    # Allelic problems identified by BLASTAlleleFinder
    IDENTITY_PROBLEMS = (32, "Allele with lower identity", None)
    SEQ_PROBLEMS = (64, "Allele sequence with non ATGC", None)
    CODON_PROBLEMS = (
        128,
        "Allele missing start or/and stop codons, or having an internal stop codon",
        None,
    )
    # Locus problem identified by BLASTAlleleFinder
    REPEAT_PROBLEMS = (256, "Locus with paralogous", None)
    # Depth problems
    COVERAGE_PROBLEMS = (
        512,
        "Allele with positions without [any nucleotide with] minimum depth",
        "CP",
    )
    BIAS_PROBLEMS = (
        1024,
        "Allele with positions with [covered nucleotide with] strand bias",
        "BP",
    )
    DOUBLE_CALLS = (2048, "Allele with positions with double possible base calls", "DC")
    MORE_CALLS = (
        4096,
        "Allele with positions with more than two possible base calls",
        "MC",
    )
    MISMATCH_CALLS = (
        8192,
        "Allele with positions with most frequent nucleotide mismatching assembly base called",
        "MM",
    )

    def __new__(cls, value: int, description: str, tag: str | None) -> Self:
        """
        Using __new__ instead of __init__ so that _value_to_member_map_
        has only the values (int) as keys and not the tuple
        """
        obj = object.__new__(cls)
        obj._value_ = value
        obj.description = description
        obj.tag = tag
        return obj

    def __init__(self, value: int, description: str, tag: str | None) -> None:
        """
        To help with auto-completion
        """
        self._value_ = value
        self.description = description
        self.tag = tag

    def __str__(self) -> str:
        return self.name

    @classmethod
    def get_flag(cls, value: int) -> Self:
        """
        Get a Flag using only the value
        """
        return cls._value2member_map_.get(value, None)


class CombinedFlag:
    def __init__(self, value: int = 0):
        self._value = value
        if value:
            self._check_flags()

    @property
    def value(self) -> int:
        self._check_flags()
        return self._value

    def set_flags(self, *flags: Flag, check_previous_assignment: bool = False):
        for flag in flags:
            if self.is_flag_set(flag):
                if check_previous_assignment:
                    raise ValueError(f"Flag {flag.name} is already set.")
            else:
                self._value ^= flag.value

    def unset_flags(self, *flags: Flag, check_previous_assignment: bool = False):
        for flag in flags:
            if not self.is_flag_set(flag):
                if check_previous_assignment:
                    raise ValueError(f"Flag {flag.name} is already unset.")
            else:
                self._value ^= flag.value

    def _check_flags(self):
        # Group
        flags_list = (Flag.CORE, Flag.ACCESSORY)
        self._check_flags_forbidden_together(flags_list)

        # Called
        flags_list = (
            Flag.CALLED,
            Flag.NOT_CALLED,
            Flag.NOT_FOUND,
        )
        self._check_flags_forbidden_together(flags_list)

        # Sequence problems
        flags_list = (
            Flag.SEQ_PROBLEMS,
            Flag.COVERAGE_PROBLEMS,
            Flag.BIAS_PROBLEMS,
            Flag.DOUBLE_CALLS,
            Flag.MORE_CALLS,
            Flag.MISMATCH_CALLS,
        )
        if any(self.is_flag_set(flag) for flag in flags_list) and (
            self.is_flag_set(Flag.CALLED) or not self.is_flag_set(Flag.NOT_CALLED)
        ):
            raise ValueError(
                f'{",".join(flag.name for flag in flags_list)} flags must be set with {Flag.NOT_CALLED} and not with {Flag.CALLED}.'
            )

    def _check_flags_mandatory_together(self, flags: list[Flag]):
        if not all(self.is_flag_set(flag) for flag in flags):
            raise ValueError(
                f'{",".join(flag.name for flag in flags)} flags must be set together.'
            )

    def _check_flags_forbidden_together(self, flags: list[Flag]):
        if sum(self.is_flag_set(flag) for flag in flags) > 1:
            raise ValueError(
                f'{",".join(flag.name for flag in flags)} flags cannot be set together.'
            )

    def is_flag_set(self, flag: Flag) -> bool:
        return self._value & flag.value != 0
