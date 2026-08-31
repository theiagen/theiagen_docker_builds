from dataclasses import asdict
from enum import StrEnum
from pathlib import Path

from ngs_pipeline_lib.tools.quality_control.assembly import SimpleFastaParser

from src.models import FlaggedAlleleCall


class CsvFields(StrEnum):
    SAMPLE_HEADER = "Sample_ID"
    NOT_CALLED_CHAR = "?"
    MULTIPLE_HITS_CONCAT_CHAR = "|"


class LocusGroups(StrEnum):
    CORE = "core"
    ACCESSORY = "accessory"


def parse_reference_loci(lines: list[str]) -> tuple[set[str], set[str]]:
    """
    Parses the <GENUS>_loci.csv file:
    ID  allele_type
    SALM_1  accessory
    SALM_2  accessory
    ...

    allele_type is either "accessory" or "core"
    returns a list for each with all the corresponding locus ID
    """
    core_ids = set()
    accessory_ids = set()

    for line in lines:
        id_, allele_type = line.strip().split("\t")
        if allele_type == LocusGroups.ACCESSORY:
            accessory_ids.add(id_)
        elif allele_type == LocusGroups.CORE:
            core_ids.add(id_)
        else:
            raise ValueError(f"Got an unexpected allele type: {allele_type}")

    return core_ids, accessory_ids


def calls_to_csv(
    sample_id: str, calls: dict[str, list[FlaggedAlleleCall]], loci: set[str]
) -> list[str]:
    rows = []
    sorted_loci = sorted(loci)
    rows.append([CsvFields.SAMPLE_HEADER] + sorted_loci)
    rows.append(
        [sample_id]
        + [
            CsvFields.MULTIPLE_HITS_CONCAT_CHAR.join(
                map(str, (allele.id for allele in calls[locus]))
            )
            if locus in calls
            else CsvFields.NOT_CALLED_CHAR
            for locus in sorted_loci
        ]
    )

    return rows


def calls_to_dict(
    sample_id: str, calls: dict[str, list[FlaggedAlleleCall]], loci: set[str]
) -> dict[str, str | dict[str, int]]:
    sorted_loci = sorted(loci)
    return {
        "sample_id": sample_id,
        "values": {
            locus: calls[locus][
                0
            ].id  # in case of multiple alleles, only report the first allele
            for locus in sorted_loci
            if calls.get(
                locus
            )  # ignore both missing locus, or locus with empty allele list
        },
    }


def calls_to_json(
    sample_id: str,
    calls: dict[str, list[FlaggedAlleleCall]],
    loci: dict[str, set[str]] = {},
) -> dict[str, str | dict[str, list[str]] | dict[str, list[FlaggedAlleleCall]]]:
    return {
        "sample_id": sample_id,
        "calls": {key: list(map(asdict, value)) for key, value in calls.items()},
        "groups": {key: list(value) for key, value in loci.items()},
    }


def get_contigs_seq(assembly: Path | str) -> dict[str, str]:
    contigs_seq: dict[str, str] = {}
    with open(assembly, encoding="utf-8") as reader:
        for contig, seq in SimpleFastaParser(reader):
            if contig in contigs_seq:
                raise ValueError(f"Repeated contig ID found: {contig}")
            contigs_seq[contig] = seq

    return contigs_seq
