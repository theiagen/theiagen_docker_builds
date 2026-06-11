from ngs_pipeline_lib.tools.tools import hash_sequence

from src.tools import format_calls


def generate_fake_core_allele_ids() -> list[str]:
    """
    Generate fake allele_ids for fake locus
    LMO_0 -> LMO_499
    """

    return format_calls(
        sample_id="sample_id",
        allele_ids={f"LMO_{i}": hash_sequence(str(i)) for i in range(500)},
    )


def generate_fake_accessory_allele_ids() -> list[str]:
    """
    Generate fake allele_ids for fake locus
    LMO_500 -> LMO_1500
    """
    return format_calls(
        sample_id="sample_id",
        allele_ids={f"LMO_{i}": hash_sequence(str(i)) for i in range(500, 1500)},
    )


def generate_fake_core_pcr_allele_ids() -> list[str]:
    """
    Generate fake allele_ids for fake locus
    LMO_2000 -> LMO_2499
    """
    return format_calls(
        sample_id="sample_id",
        allele_ids={f"LMO_{i}": hash_sequence(str(i)) for i in range(2000, 2500)},
    )


def generate_fake_accessory_pcr_allele_ids() -> list[str]:
    """
    Generate fake allele_ids for fake locus
    LMO_2500 -> LMO_3500
    """
    return format_calls(
        sample_id="sample_id",
        allele_ids={f"LMO_{i}": hash_sequence(str(i)) for i in range(2500, 3500)},
    )


def generate_fake_result() -> dict:
    """
    Generate fake wgmlst profile
    {LMO_0: 1 -> LMO_499: 1}
    """

    return {"values": {f"LMO_{i}": "1" for i in range(500)}}
