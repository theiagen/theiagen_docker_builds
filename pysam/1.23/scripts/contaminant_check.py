import sys
import json
import logging
import argparse
from collections import defaultdict

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def write_json(filename, data):
    """Write a JSON file compatible with WDL"""
    with open(filename, "w") as f:
        if data:
            json.dump(data, f, indent=4)
        else:
            # spoof Cromwell (Terra WDL)
            f.write('{"": 0}')


def apply_thresholds(variable_by_sequence, min_value):
    """Apply thresholds for a variable while excluding empty values as missing"""
    passing_sequences_prep = set(
        [seq for seq, val in variable_by_sequence.items() if val >= min_value]
    )
    failing_sequences = set(
        [seq for seq, val in variable_by_sequence.items() if not val > 0]
    )
    passing_sequences = passing_sequences_prep.difference(failing_sequences)
    return failing_sequences, passing_sequences


def compile_failures(
    passing_sequences, expected_recovered_sequences, variable_by_sequence, var_name
):
    """Compile failing sequences based on those that are expected v. unexpected"""
    # sequences that passed coverage threshold and were expected
    hit_sequences = passing_sequences.intersection(expected_recovered_sequences)
    if hit_sequences:
        logger.debug(f"passing {var_name}: {sorted(hit_sequences)}")
    else:
        logger.warning(f"no sequences passed {var_name} threshold")
    # sequences that failed variable threshold and were expected
    missing_sequences = expected_recovered_sequences.difference(hit_sequences)
    if missing_sequences:
        logger.warning(f"failing {var_name}: {sorted(missing_sequences)}")
    # sequences that passed variable threshold but were not expected
    extra_sequences = passing_sequences.difference(expected_recovered_sequences)
    unexpected_sequences = {
        seq: variable_by_sequence[seq] for seq in sorted(extra_sequences)
    }
    write_json(f"UNEXPECTED_SEQ2{var_name.upper()}.json", unexpected_sequences)
    return missing_sequences, extra_sequences


def annotate_failures(
    seq2fail, variable_missing, failing_sequences, variable_by_sequence, variable
):
    """Annotate failure type for the status output"""
    for seq in variable_missing:
        if seq in failing_sequences or seq not in variable_by_sequence:
            seq2fail[seq].append("not detected")
        else:
            seq2fail[seq].append(
                f"insufficient {variable} ({variable_by_sequence[seq]})"
            )
    return seq2fail


def write_status(
    expected_recovered_sequences,
    minimum_expected_sequences,
    unexpected_sequences,
    maximum_unexpected_sequences,
    seq2fail,
):
    # populate a status string
    with open("STATUS", "w") as f:
        # check if a pass/fail threshold was infringed
        if (
            len(expected_recovered_sequences) < minimum_expected_sequences
            or len(unexpected_sequences) > maximum_unexpected_sequences
        ):
            status_string = "FAIL: "
            # too few expected sequences recovered
            if len(expected_recovered_sequences) < minimum_expected_sequences:
                for seq, fail_reasons in sorted(seq2fail.items(), key=lambda x: x[0]):
                    status_string += f"{seq} - {', '.join(fail_reasons)}; "
            # too many unexpected sequences recovered
            if len(unexpected_sequences) > maximum_unexpected_sequences:
                for seq in unexpected_sequences:
                    status_string += f"{seq} - extra sequence; "
            status_string = status_string.strip("; ")
            f.write(status_string)
        else:
            f.write("PASS")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="task_contaminant_check.py dependency script (github.com/theiagen/public_health_bioinformatics)"
    )
    parser.add_argument("--expected_sequences", required=True)
    parser.add_argument("--minimum_expected_sequences")
    parser.add_argument("--maximum_unexpected_sequences", default=0)
    parser.add_argument("--coverage_by_sequence_json", required=True)
    parser.add_argument("--depth_by_sequence_json", required=True)
    parser.add_argument("--reads_by_sequence_json", required=True)
    parser.add_argument("--contaminant_fasta", required=True)
    parser.add_argument("--minimum_percent_coverage", required=True)
    parser.add_argument("--minimum_depth", required=True)
    parser.add_argument("--minimum_reads_mapped", required=True)
    args = parser.parse_args()

    # convert comma-separated string of expected sequences into a set
    expected_sequences = set(
        [
            seq.strip()
            for seq in args.expected_sequences.strip('"').strip("'").split(",")
            if seq.strip()
        ]
    )
    # set default to all expected_sequences
    if args.minimum_expected_sequences is not None:
        minimum_expected_sequences = args.minimum_expected_sequences
        if minimum_expected_sequences > len(expected_sequences):
            logger.error(
                f"minimum expected sequences ({minimum_expected_sequences}) exceeds number of expected sequences ({len(expected_sequences)}); setting minimum expected sequences to {len(expected_sequences)}"
            )
            minimum_expected_sequences = len(expected_sequences)
    else:
        minimum_expected_sequences = len(expected_sequences)
    logger.debug(f"expecting minimum {minimum_expected_sequences} sequences")

    # read in coverage and depth by sequence
    with open(args.coverage_by_sequence_json) as f:
        coverage_by_sequence = json.load(f)
    with open(args.depth_by_sequence_json) as f:
        depth_by_sequence = json.load(f)
    with open(args.reads_by_sequence_json) as f:
        reads_by_sequence = json.load(f)
    # obtain all sequences in the reference
    reference_sequences = set()
    with open(args.contaminant_fasta) as f:
        for line in f:
            if line.startswith(">"):
                # acquire the sequence from the header by removing descriptions and the ">"
                seq = line.split()[0][1:]
                reference_sequences.add(seq)

    # check if any expected sequences are present above the specified thresholds
    failing_sequences_coverage, passing_sequences_coverage = apply_thresholds(
        coverage_by_sequence, args.mininimum_percent_coverage
    )
    failing_sequences_depth, passing_sequences_depth = apply_thresholds(
        depth_by_sequence, args.minimum_depth
    )
    failing_sequences_reads, passing_sequences_reads = apply_thresholds(
        reads_by_sequence, args.mininimum_reads_mapped
    )

    failing_sequences = failing_sequences_coverage.union(failing_sequences_depth).union(
        failing_sequences_reads
    )
    passing_sequences = passing_sequences_coverage.union(passing_sequences_depth).union(
        passing_sequences_reads
    )

    # sequences that were desired to be identified, but not recovered
    expected_unrecovered_sequences = expected_sequences.difference(reference_sequences)
    # sequences that were expected and recovered
    expected_recovered_sequences = expected_sequences.difference(
        expected_unrecovered_sequences
    )

    # write outputs for recovered expected sequences
    expected_sequences_coverage = {}
    expected_sequences_depth = {}
    expected_sequences_reads = {}
    for seq in sorted(expected_recovered_sequences):
        expected_sequences_coverage[seq] = coverage_by_sequence[seq]
        expected_sequences_depth[seq] = depth_by_sequence[seq]
        expected_sequences_reads[seq] = reads_by_sequence[seq]
    write_json("EXPECTED_SEQ2COVERAGE.json", expected_sequences_coverage)
    write_json("EXPECTED_SEQ2DEPTH.json", expected_sequences_depth)
    write_json("EXPECTED_SEQ2READS.json", expected_sequences_reads)

    # check results and write outputs for unexpected sequences
    logger.debug(f"expected sequences: {sorted(expected_sequences)}")
    coverage_missing, coverage_extra = compile_failures(
        passing_sequences_coverage,
        expected_recovered_sequences,
        coverage_by_sequence,
        "coverage",
    )
    depth_missing, depth_extra = compile_failures(
        passing_sequences_depth,
        expected_recovered_sequences,
        depth_by_sequence,
        "depth",
    )
    reads_missing, reads_extra = compile_failures(
        passing_sequences_reads,
        expected_recovered_sequences,
        reads_by_sequence,
        "reads",
    )

    # total unexpected sequences that passed thresholds
    unexpected_sequences = sorted(coverage_extra.union(depth_extra).union(reads_extra))
    if unexpected_sequences:
        logger.warning(f"extraneous sequences detected: {sorted(unexpected_sequences)}")

    # annotate failing sequences
    seq2fail = defaultdict(list)
    seq2fail = annotate_failures(
        seq2fail, coverage_missing, failing_sequences, coverage_by_sequence, "coverage"
    )
    seq2fail = annotate_failures(
        seq2fail, depth_missing, failing_sequences, depth_by_sequence, "depth"
    )
    seq2fail = annotate_failures(
        seq2fail, reads_missing, failing_sequences, reads_by_sequence, "reads mapped"
    )
    # these sequences are missing from the reference because they were expected
    # but not detected/accounted for in the reference FASTA
    for seq in expected_unrecovered_sequences:
        seq2fail[seq].append("missing from reference")

    write_status(
        expected_recovered_sequences,
        minimum_expected_sequences,
        unexpected_sequences,
        args.maximum_unexpected_sequences,
        seq2fail,
    )

    sys.exit(0)
