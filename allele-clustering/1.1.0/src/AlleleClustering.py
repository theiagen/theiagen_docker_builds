#!/usr/bin/env python3

import argparse
import csv
import json

from distance_matrix import DistanceMatrix
from typing import Any


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments

    Returns:
        argparse.Namespace: Parsed command-line arguments containing the
        input allele profile file, clsutering algorithm, distance
        metric, and output file prefix
    """
    parser = argparse.ArgumentParser(
        description="Generate a NWK from allele profiles in NDJSON format."
    )

    parser.add_argument(
        "allele_profiles",
        help="The input NDJSON file containing newline-concatenated allele profiles.",
    )
    parser.add_argument(
        "-a",
        "--algorithm",
        choices=[
            "upgma",
            "single_linkage",
            "complete_linkage",
            "neighbor_joining",
            "minimum_spanning",
        ],
        required=True,
        help="The algorithm to use. `minimum_spanning` is an implemntation intended to reproduce the behavior of the BioNumerics algorithm, while `neighbor_joining` uses a dynamic neighbor joining algorithm. `upgma`, `single_linkage`, and `complete_linkage` are also available options",
    )
    parser.add_argument(
        "-d",
        "--distance",
        choices=["absolute_allele_differences", "normalized_allele_differences"],
        required=True,
        help="Distance metric to use. `absolute_allele_differences` uses a raw count of the positions at which a pair of samples differ (missing data is not counted as a difference), while `normalized_allele_differences` uses the number of positions at which two samples differ divided by the number of possitions for which both have data",
    )

    parser.add_argument(
        "-o",
        "--output",
        help="The prefix for the generated NWK file",
        default="allele_clustering",
    )

    return parser.parse_args()


def parse_json(ndjson_file) -> dict[str, Any]:
    """Convert an NDJSON allele profile file into an allele matrix.

    Each line of the input file must be a JSON object containing a
    ``sample_id`` field and a ``values`` dictionary mapping locus names
    to allele identifiers. All loci observed across samples are collected,
    sorted alphabetically, and used to construct a matrix where each row
    corresponds to a sample and missing loci are represented by an empty
    string.

    Args:
        ndjson_file: Path to the input NDJSON file.

    Returns:
        dict[str, Any]: Dictionary containing the allele matrix. The
            ``"Headers"`` key stores the ordered list of loci, while each
            remaining key is a sample ID whose value is a list of allele
            values aligned with the headers.
    """
    rows = {}
    all_loci = set()

    with open(ndjson_file, "r") as infile:
        for line in infile:
            sample_record = json.loads(line)

            sample_id = sample_record["sample_id"]
            hashes = sample_record["values"]

            rows[sample_id] = hashes
            all_loci.update(hashes.keys())

    # not sure if sorting is needed but it doesn't hurt
    sorted_loci = sorted(all_loci)

    allele_matrix = {"Headers": sorted_loci}

    for sample_id, hashes in rows.items():
        allele_matrix[sample_id] = [hashes.get(locus, "") for locus in sorted_loci]

    return allele_matrix


if __name__ == "__main__":
    args = parse_arguments()

    # generate tree
    allele_profile_matrix = parse_json(args.allele_profiles)
    tree = DistanceMatrix.from_json(
        allele_profile_matrix, algorithm=args.algorithm, distance=args.distance
    )

    # write nwk tree to file
    output_name = args.output + ".nwk"
    with open(output_name, "w") as output:
        output.write(tree.tree)

    # write allele distance matrix to a csv file
    output_name = args.output + ".matrix.csv"
    with open(output_name, "w") as output:
        output_matrix = csv.writer(output)
        output_matrix.writerows(tree.dmat_as_list())
