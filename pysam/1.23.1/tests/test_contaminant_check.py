from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "contaminant_check.py"
SPEC = spec_from_file_location("contaminant_check", SCRIPT_PATH)
contaminant_check = module_from_spec(SPEC)
SPEC.loader.exec_module(contaminant_check)


def _run_contaminant_check(
    tmp_path,
    coverage_json,
    depth_json,
    reads_json,
    fasta_text,
    expected_sequences="seq1,seq2",
    extra_args=None,
    min_percent_coverage=90,
    min_depth=5,
    min_reads_mapped=5
):
    coverage = tmp_path / "coverage.json"
    depth = tmp_path / "depth.json"
    reads = tmp_path / "reads.json"
    fasta = tmp_path / "contaminants.fasta"

    coverage.write_text(coverage_json)
    depth.write_text(depth_json)
    reads.write_text(reads_json)
    fasta.write_text(fasta_text)

    cmd = [
        sys.executable,
        str(SCRIPT_PATH),
        "--expected_sequences",
        expected_sequences,
        "--coverage_by_sequence_json",
        str(coverage),
        "--depth_by_sequence_json",
        str(depth),
        "--reads_by_sequence_json",
        str(reads),
        "--contaminant_fasta",
        str(fasta),
        "--minimum_percent_coverage",
        str(min_percent_coverage),
        "--minimum_depth",
        str(min_depth),
        "--minimum_reads_mapped",
        str(min_reads_mapped),
    ]
    if extra_args:
        cmd.extend(extra_args)

    subprocess.run(cmd, cwd=tmp_path, check=True)


def test_threshold_gate_requires_all_three_metrics(tmp_path):
    _run_contaminant_check(
        tmp_path,
        '{"seq1": 100.0, "seq2": 100.0}',
        '{"seq1": 20, "seq2": 0}',
        '{"seq1": 50, "seq2": 10}',
        ">seq1\nACTG\n>seq2\nACTG\n",
    )

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text.startswith("FAIL:")
    assert "seq2" in status_text


def test_passing_sample_with_mixed_expected_and_allowed_unexpected(tmp_path):
    _run_contaminant_check(
        tmp_path,
        '{"seqA": 95.0, "seqB": 90.0, "seqC": 0, "seqX": 91.0}',
        '{"seqA": 10, "seqB": 5, "seqC": 0, "seqX": 7}',
        '{"seqA": 30, "seqB": 5, "seqC": 0, "seqX": 9}',
        ">seqA\nACTG\n>seqB\nACTG\n>seqX\nACTG\n",
        expected_sequences="seqA,seqB,seqC",
        extra_args=[
            "--minimum_expected_sequences",
            "2",
            "--maximum_unexpected_sequences",
            "1",
        ],
    )

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text == "PASS"


def test_fails_when_minimum_expected_not_met(tmp_path):
    _run_contaminant_check(
        tmp_path,
        '{"seq1": 100.0, "seq2": 100.0, "seq3": 89.0}',
        '{"seq1": 20, "seq2": 8, "seq3": 2}',
        '{"seq1": 50, "seq2": 12, "seq3": 2}',
        ">seq1\nACTG\n>seq2\nACTG\n>seq3\nACTG\n",
        expected_sequences="seq1,seq2,seq3",
        extra_args=["--minimum_expected_sequences", "3"],
    )

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text.startswith("FAIL:")
    assert "seq3" in status_text


def test_fails_when_maximum_unexpected_exceeded(tmp_path):
    _run_contaminant_check(
        tmp_path,
        '{"seq1": 100.0, "seq2": 100.0, "seqX": 95.0, "seqY": 99.0}',
        '{"seq1": 20, "seq2": 20, "seqX": 10, "seqY": 11}',
        '{"seq1": 50, "seq2": 50, "seqX": 9, "seqY": 12}',
        ">seq1\nACTG\n>seq2\nACTG\n>seqX\nACTG\n>seqY\nACTG\n",
        expected_sequences="seq1,seq2",
        extra_args=["--maximum_unexpected_sequences", "1"],
    )

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text.startswith("FAIL:")
    assert "extra sequence" in status_text


def test_fails_when_expected_sequence_missing_from_reference(tmp_path):
    _run_contaminant_check(
        tmp_path,
        '{"seq1": 100.0, "seq2": 100.0}',
        '{"seq1": 20, "seq2": 20}',
        '{"seq1": 50, "seq2": 50}',
        ">seq1\nACTG\n",
        expected_sequences="seq1,seq2",
    )

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text.startswith("FAIL:")
    assert "seq2 - missing from reference" in status_text


def test_passes_when_unexpected_count_equals_maximum(tmp_path):
    _run_contaminant_check(
        tmp_path,
        '{"seq1": 100.0, "seq2": 100.0, "seqX": 95.0}',
        '{"seq1": 20, "seq2": 20, "seqX": 10}',
        '{"seq1": 50, "seq2": 50, "seqX": 9}',
        ">seq1\nACTG\n>seq2\nACTG\n>seqX\nACTG\n",
        expected_sequences="seq1,seq2",
        extra_args=["--maximum_unexpected_sequences", "1"],
    )

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text == "PASS"
