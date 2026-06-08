from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "contaminant_check.py"
SPEC = spec_from_file_location("contaminant_check", SCRIPT_PATH)
contaminant_check = module_from_spec(SPEC)
SPEC.loader.exec_module(contaminant_check)


def test_threshold_gate_requires_all_three_metrics(tmp_path):
    coverage = tmp_path / "coverage.json"
    depth = tmp_path / "depth.json"
    reads = tmp_path / "reads.json"
    fasta = tmp_path / "contaminants.fasta"

    coverage.write_text('{"seq1": 100.0, "seq2": 100.0}')
    depth.write_text('{"seq1": 20, "seq2": 0}')
    reads.write_text('{"seq1": 50, "seq2": 10}')
    fasta.write_text(">seq1\nACTG\n>seq2\nACTG\n")

    cmd = [
        sys.executable,
        str(SCRIPT_PATH),
        "--expected_sequences",
        "seq1,seq2",
        "--coverage_by_sequence_json",
        str(coverage),
        "--depth_by_sequence_json",
        str(depth),
        "--reads_by_sequence_json",
        str(reads),
        "--contaminant_fasta",
        str(fasta),
        "--minimum_percent_coverage",
        "90",
        "--minimum_depth",
        "5",
        "--minimum_reads_mapped",
        "5",
    ]

    subprocess.run(cmd, cwd=tmp_path, check=True)

    status_text = (tmp_path / "STATUS").read_text().strip()
    assert status_text.startswith("FAIL:")
    assert "seq2" in status_text
