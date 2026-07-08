#!/usr/bin/env python3

from distance_matrix import DistanceMatrix

eg_profile = {
    "Headers": ["1", "HEY", "3", "4", "5"],
    "A": ["1", "1", "1", "1", "1"],
    "B": ["1", "1", "1", "1", "2"],
    "CAT": ["2", "2", "2", "3", "3"],
    "D": ["2", "", "3", "1", "3"],
    "E": ["4", "4", "4", "4", "4"],
}

dmat = DistanceMatrix.from_json(
    eg_profile, algorithm="neighbor_joining", distance="absolute_allele_differences"
)
print(dmat.tree)
print(dmat.dmat_as_list())
