# PulseNet 2.0 Allele Clustering container

**Main tool**: [PulseNet 2.0 AlleleClustering](https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main)

**Code repository**: <https://github.com/ncezid-biome/pulsenet2.0-trees/tree/main>

This directory copies the following files from the aforementioned GitHub repository:
- [`distance_matrix.py`](https://github.com/ncezid-biome/pulsenet2.0-trees/blob/main/distance_matrix.py)
- [`example.py`](https://github.com/ncezid-biome/pulsenet2.0-trees/blob/main/example.py)
- [`tree_model.py`](https://github.com/ncezid-biome/pulsenet2.0-trees/blob/main/tree_model.py)

**Basic information on how to use this tool**:
- executable: `python3 AlleleClustering.py`
- help: `python3 AlleleClustering.py --help`
- description: the PulseNet 2.0 AlleleClustering algorithm

**Additional information**:

- contains numpy==2.4.4

**Full documentation**: [See the README.md here](https://github.com/ncezid-biome/pulsenet2.0-trees/blob/main/README.md)

## Example Usage

This example shows how to ensure that the Python code `AlleleClustering.py` matches the native implementation of the algorithm shown in `example.py`

First, establish the "truth" by running the following command (assuming you are in the 1.0.0 directory):

```
$ python3 src/example.py
E:2.625,((A:0.5,B:0.5):1.75,(D:0.75,CAT:1.25):1.25):0.125);
```

Then, to confirm proper functioning, run the following command from the 1.0.0 directory (or alter paths):

```
$ python3 src/AlleleClustering.py ../tests/test.ndjson -a neighbor_joining -d absolute_allele_differences && cat allele_clustering.nwk && echo ""
(E:2.625,((A:0.5,B:0.5):1.75,(D:0.75,CAT:1.25):1.25):0.125);
```

## Changelog

- v1.0.0 copies the algorithm files from the home repository and creates the `AlleleClustering.py` interface to transform NDJSONs into a format useable for the algorithm.
- v1.1.0 alters `AlleleClustering.py` by adding docstrings and returning the results from the `DistanceMatrix.from_json` result's `dmat_as_list()` function to a CSV file.
- v1.2.0 alters `AlleleClustering.py` by adding column and row names to the distance matrix and enables this file to be generated when the algorithm is "minimum_spanning"
