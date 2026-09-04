from pytest import mark, raises

from src.sam_flags import CombinedFlag, Flag


def test_flag():
    assert Flag.CORE.value == 1
    assert Flag.CORE.name == "CORE"
    assert str(Flag.CORE) == "CORE"
    assert Flag.CORE.description == "Core locus"

    assert len(Flag._value2member_map_.keys()) == 14
    assert Flag.get_flag(1) == Flag.CORE
    assert Flag.get_flag(4096) == Flag.MORE_CALLS


def test_combined_flag_valid_init():
    assert CombinedFlag(value=10)


@mark.parametrize(
    "value",
    [
        3,  # CORE + ACCESSORY -> forbidden together
        12,  # CALLED + NOT_CALLED -> forbidden_together
        20,  # CALLED + NOT_FOUND -> forbidden_together
        24,  # NOT_CALLED + NOT_FOUND -> forbidden_together
        64,  # SEQ_PROBLEMS -> must have NOT_CALLED
        68,  # SEQ_PROBLEMS + CALLED -> can't be CALLED
        512,  # COVERAGE_PROBLEMS -> must have NOT_CALLED
        516,  # COVERAGE_PROBLEMS + CALLED -> can't be CALLED
        1024,  # BIAIS_PROBLEMS -> must have NOT_CALLED
        1028,  # BIAIS_PROBLEMS + CALLED -> can't be CALLED
        2048,  # DOUBLE_CALLS -> must have NOT_CALLED
        2052,  # DOUBLE_CALLS + CALLED -> can't be CALLED
        4096,  # MORE_CALLS -> must have NOT_CALLED
        4100,  # MORE_CALLS + CALLED -> can't be CALLED
        8192,  # MISMATCH_CALLS -> must have NOT_CALLED
        8196,  # MISMATCH_CALLS + CALLED -> can't be CALLED
    ],
)
def test_combined_flag_invalid_init(value: int):
    with raises(ValueError):
        CombinedFlag(value=value)


def test_combined_flag_valid_set_flags():
    combined_flag = CombinedFlag()
    combined_flag.set_flags(Flag.CORE, Flag.CALLED)
    assert combined_flag.value == 5
    # already set does not change
    combined_flag.set_flags(Flag.CORE)
    assert combined_flag.value == 5


def test_combined_flag_invalid_set_flags():
    combined_flag = CombinedFlag(value=1)
    with raises(ValueError):
        combined_flag.set_flags(Flag.CORE, check_previous_assignment=True)


def test_combined_flag_valid_unset_flags():
    combined_flag = CombinedFlag(value=5)
    combined_flag.unset_flags(Flag.CORE)
    assert combined_flag.value == 4
    # already unset does not change
    combined_flag.unset_flags(Flag.CORE)
    assert combined_flag.value == 4


def test_combined_flag_invalid_unset_flags():
    combined_flag = CombinedFlag(value=2)
    with raises(ValueError):
        combined_flag.unset_flags(Flag.CORE, check_previous_assignment=True)


def test_combined_flag_valid_check_mandatory_together():
    combined_flag = CombinedFlag(value=5)
    combined_flag._check_flags_mandatory_together([Flag.CORE, Flag.CALLED])


def test_combined_flag_invalid_check_mandatory_together():
    combined_flag = CombinedFlag(value=1)
    with raises(ValueError):
        combined_flag._check_flags_mandatory_together([Flag.CORE, Flag.CALLED])


def test_combined_flag_valid_check_forbidden_together():
    combined_flag = CombinedFlag(value=6)
    combined_flag._check_flags_forbidden_together([Flag.CORE, Flag.CALLED])


def test_combined_flag_invalid_check_forbidden_together():
    combined_flag = CombinedFlag(value=5)
    with raises(ValueError):
        combined_flag._check_flags_forbidden_together([Flag.CORE, Flag.CALLED])
