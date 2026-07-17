from collections import defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import re

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "gene_coverage.py"
SPEC = spec_from_file_location("gene_coverage", SCRIPT_PATH)
gene_coverage = module_from_spec(SPEC)
SPEC.loader.exec_module(gene_coverage)


class MockBam:
    def __init__(self, references, contig_lengths, default_depth=5):
        self.references = tuple(references)
        self._lengths = dict(contig_lengths)
        self._default_depth = default_depth

    def get_reference_length(self, contig):
        return self._lengths[contig]

    def count_coverage(self, contig, start, end, quality_threshold=0):
        span = end - start
        return (
            [self._default_depth] * span,
            [0] * span,
            [0] * span,
            [0] * span,
        )


def _write_mock_gbff(path, contig="contig1", gene="geneA", start=10, end=20):
    record = SeqRecord(Seq("ATCG" * 50), id=contig, name=contig, description="")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(
        SeqFeature(
            FeatureLocation(start, end),
            type="CDS",
            qualifiers={"product": [gene]},
        )
    )
    with open(path, "w") as handle:
        SeqIO.write(record, handle, "genbank")


def test_parse_gbff_extracts_expected_coordinates(tmp_path):
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 100}, default_depth=5
    )
    gbff = tmp_path / "mock.gbff"
    _write_mock_gbff(gbff, contig="contig1", gene="geneA", start=10, end=20)

    contig2query2coords = defaultdict(lambda: defaultdict(list))
    parsed = gene_coverage.parse_gbff(
        set(mock_bam.references),
        str(gbff),
        {"geneA"},
        "CDS",
        "product",
        gene_coverage.exact_check,
        contig2query2coords,
        False
    )

    assert parsed["contig1"]["geneA"] == [[10, 20]]


def test_bed_and_gbff_coordinates_agree_for_same_gene(tmp_path):
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 100}, default_depth=5
    )

    gbff = tmp_path / "mock.gbff"
    _write_mock_gbff(gbff, contig="contig1", gene="geneA", start=10, end=20)

    bed = tmp_path / "mock.bed"
    bed.write_text("contig1\t10\t20\tgeneA\n")

    from_gbff = gene_coverage.parse_gbff(
        set(mock_bam.references),
        str(gbff),
        {"geneA"},
        "CDS",
        "product",
        gene_coverage.exact_check,
        defaultdict(lambda: defaultdict(list)),
        False
    )
    from_bed = gene_coverage.parse_bed(
        set(mock_bam.references),
        str(bed),
        {"geneA"},
        gene_coverage.exact_check,
        defaultdict(lambda: defaultdict(list)),
        False
    )

    gbff_coords = [tuple(x) for x in from_gbff["contig1"]["geneA"]]
    bed_coords = [tuple(x) for x in from_bed["contig1"]["geneA"]]
    assert gbff_coords == bed_coords == [(10, 20)]


def _write_mock_vcf(path, contig="contig1", contig_length=100, variants=None):
    """Write a minimal VCF; variants is a list of 1-based POS integers"""
    variants = variants or []
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={contig},length={contig_length}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    lines.extend(f"{contig}\t{pos}\t.\tA\tT\t.\t.\t." for pos in variants)
    path.write_text("\n".join(lines) + "\n")


def test_extract_vcf_genes_indexing_is_zero_based_half_open(tmp_path):
    import pysam

    # geneA occupies 0-based half-open [10, 20) == 1-based positions 11..20
    vcf_in = tmp_path / "in.vcf"
    _write_mock_vcf(
        vcf_in, contig="contig1", contig_length=100, variants=[10, 11, 20, 21]
    )
    output = tmp_path / "out.vcf"

    contig2query2coords = {"contig1": {"geneA": [(10, 20)]}}
    written = gene_coverage.extract_vcf_genes(
        str(vcf_in), contig2query2coords, str(output)
    )

    with pysam.VariantFile(str(output)) as handle:
        assert "GENE" in handle.header.info
        kept = [(rec.pos, tuple(rec.info["GENE"])) for rec in handle]

    # POS 10 (base index 9) sits before the range; POS 21 (base index 20) sits past it.
    # POS 11 (first base) and POS 20 (last base) fall within [10, 20).
    assert written == 2
    assert kept == [(11, ("geneA",)), (20, ("geneA",))]


def test_quantify_gene_coverage_known_depth_and_breadth():
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 100}, default_depth=5
    )
    contig2query2coords = {"contig1": {"geneA": [(10, 20)]}}

    depth_dict, coverage_dict = gene_coverage.quantify_gene_coverage(
        mock_bam,
        contig2query2coords,
        min_depth=1,
        min_quality=0,
    )

    assert depth_dict["geneA"] == 5.0
    assert coverage_dict["geneA"] == 100.0


@pytest.mark.parametrize(
    "coords, message_fragment",
    [
        ({"contig1": {"geneA": [(10, 10)]}}, "start (10) must be < end (10)"),
        (
            {"contig1": {"geneA": [(5, 15)]}},
            "end (15) exceeds contig length (10)",
        ),
        (
            {"missing_contig": {"geneA": [(1, 5)]}},
            "not found in BAM references",
        ),
    ],
)
def test_quantify_gene_coverage_edge_guards_raise_clean_value_errors(
    coords, message_fragment
):
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 10}, default_depth=5
    )

    with pytest.raises(ValueError, match=re.escape(message_fragment)):
        gene_coverage.quantify_gene_coverage(
            mock_bam, coords, min_depth=1, min_quality=0
        )
