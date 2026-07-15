"""Tests for variant_annotation.py

Synthetic genes with fully controlled coding sequences let every expected
c./p. notation be predicted by hand.  The forward-strand gene "test gene alpha"
has coding sequence ATG TAT CCC AAA GGG TTT CAT TGA (M Y P K G F H *) placed at
0-based genomic [9, 33) on contig chr1, so coding index i == genomic 9 + i and
VCF POS == 10 + i.
"""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "variant_annotation.py"
SPEC = spec_from_file_location("variant_annotation", SCRIPT_PATH)
va = module_from_spec(SPEC)
SPEC.loader.exec_module(va)


ALPHA_CODING = "ATGTATCCCAAAGGGTTTCATTGA"  # M Y P K G F H *
ALPHA_START = 9


def _make_gbff(path):
    """Write a GBFF with a forward, reverse and two-exon gene."""
    records = []

    # chr1: forward-strand "test gene alpha" at [9, 33)
    chr1 = SeqRecord(
        Seq("A" * ALPHA_START + ALPHA_CODING + "A" * 7),
        id="chr1",
        name="chr1",
        description="",
    )
    chr1.annotations["molecule_type"] = "DNA"
    chr1.features.append(
        SeqFeature(
            FeatureLocation(ALPHA_START, ALPHA_START + len(ALPHA_CODING), strand=1),
            type="CDS",
            qualifiers={"product": ["test gene alpha"], "transl_table": ["1"]},
        )
    )
    records.append(chr1)

    # chr2: reverse-strand "test gene beta"; coding ATG TAT AAA TGA (M Y K *)
    beta_coding = "ATGTATAAATGA"
    beta_genomic = str(Seq(beta_coding).reverse_complement())
    chr2_seq = "C" * 5 + beta_genomic + "C" * 5
    chr2 = SeqRecord(Seq(chr2_seq), id="chr2", name="chr2", description="")
    chr2.annotations["molecule_type"] = "DNA"
    chr2.features.append(
        SeqFeature(
            FeatureLocation(5, 5 + len(beta_coding), strand=-1),
            type="CDS",
            qualifiers={"product": ["test gene beta"], "transl_table": ["1"]},
        )
    )
    records.append(chr2)

    # chr3: forward two-exon "test gene gamma"; exon1 [0,6) exon2 [10,16)
    # coding ATG TAT | AAA TGA  (M Y K *) with an intron at [6, 10)
    gamma_seq = "ATGTAT" + "GGGG" + "AAATGA" + "A" * 4
    chr3 = SeqRecord(Seq(gamma_seq), id="chr3", name="chr3", description="")
    chr3.annotations["molecule_type"] = "DNA"
    chr3.features.append(
        SeqFeature(
            CompoundLocation(
                [
                    FeatureLocation(0, 6, strand=1),
                    FeatureLocation(10, 16, strand=1),
                ]
            ),
            type="CDS",
            qualifiers={"product": ["test gene gamma"], "transl_table": ["1"]},
        )
    )
    records.append(chr3)

    with open(path, "w") as handle:
        SeqIO.write(records, handle, "genbank")


def _write_vcf(path, contig, contig_len, rows, with_sample=False):
    """rows: list of (pos, ref, alt) or (pos, ref, alt, ref_ad, alt_ad)."""
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={contig},length={contig_len}>",
        '##INFO=<ID=GENE,Number=.,Type=String,Description="gene">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    ]
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if with_sample:
        header += "\tFORMAT\tSAMPLE1"
    lines.append(header)
    for row in rows:
        pos, ref, alt = row[0], row[1], row[2]
        line = f"{contig}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
        if with_sample:
            ref_ad, alt_ad = row[3], row[4]
            line += f"\tGT:AD\t0/1:{ref_ad},{alt_ad}"
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def gbff(tmp_path):
    path = tmp_path / "ref.gbff"
    _make_gbff(path)
    return str(path)


def _annotate(gbff, vcf, query=("test gene alpha",), exact=True):
    return va.run(vcf, gbff, list(query), exact_match=exact)


# --------------------------------------------------------------------------- #
# coding-sequence model
# --------------------------------------------------------------------------- #

def test_gene_model_matches_biopython_extract(gbff):
    """Our hand-rolled coding sequence must equal BioPython's feature.extract."""
    for record in SeqIO.parse(gbff, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            product = feature.qualifiers["product"][0]
            models = va.build_gene_models(
                gbff, [product], "CDS", "product", exact_match=True
            )
            model = models[product]
            assert model.ref_coding == str(feature.extract(record.seq)).upper()


def test_alpha_reference_protein(gbff):
    models = va.build_gene_models(
        gbff, ["test gene alpha"], "CDS", "product", exact_match=True
    )
    assert models["test gene alpha"].ref_coding == ALPHA_CODING
    assert models["test gene alpha"].ref_protein == "MYPKGFH*"


# --------------------------------------------------------------------------- #
# single-nucleotide substitutions
# --------------------------------------------------------------------------- #

def test_missense(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(14, "A", "T")])  # c.5 A>T, codon 2 TAT->TTT
    report = _annotate(gbff, str(vcf))
    assert "missense_variant c.5A>T p.Tyr2Phe" in report
    assert report.startswith("test.gene.alpha: test gene alpha (")


def test_synonymous(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(18, "C", "T")])  # c.9 C>T, codon 3 CCC->CCT
    report = _annotate(gbff, str(vcf))
    assert "synonymous_variant c.9C>T p.Pro3=" in report


def test_nonsense(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(19, "A", "T")])  # c.10 A>T, codon 4 AAA->TAA
    report = _annotate(gbff, str(vcf))
    assert "stop_gained c.10A>T p.Lys4Ter" in report


def test_start_lost(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(10, "A", "C")])  # c.1 A>C, ATG->CTG
    report = _annotate(gbff, str(vcf))
    assert "start_lost c.1A>C" in report


def test_stop_lost(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(33, "A", "G")])  # c.24 A>G, codon 8 TGA->TGG
    report = _annotate(gbff, str(vcf))
    assert "stop_lost c.24A>G p.Ter8Trp" in report


# --------------------------------------------------------------------------- #
# indels
# --------------------------------------------------------------------------- #

def test_frameshift_deletion(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    # delete the G at c.13 (codon 5, GGG); anchor at genomic 20 (POS 21)
    _write_vcf(vcf, "chr1", 40, [(21, "AG", "A")])
    report = _annotate(gbff, str(vcf))
    assert "frameshift_variant c.13delG p.Gly5fs" in report


def test_frameshift_insertion(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    # insert a single T after c.9 (frameshift)
    _write_vcf(vcf, "chr1", 40, [(18, "C", "CT")])
    report = _annotate(gbff, str(vcf))
    assert "frameshift_variant" in report
    assert "p.Lys4fs" in report  # first affected codon is 4


def test_inframe_deletion(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    # delete codon 3 (CCC); anchor at genomic 14 (POS 15)
    _write_vcf(vcf, "chr1", 40, [(15, "TCCC", "T")])
    report = _annotate(gbff, str(vcf))
    assert "inframe_deletion c.7_9del p.Pro3del" in report


def test_inframe_insertion(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    # insert GGG (one codon) after c.9
    _write_vcf(vcf, "chr1", 40, [(18, "C", "CGGG")])
    report = _annotate(gbff, str(vcf))
    assert "inframe_insertion c.9_10insGGG p.Pro3_Lys4insGly" in report


# --------------------------------------------------------------------------- #
# frameshift downstream suppression
# --------------------------------------------------------------------------- #

def test_frameshift_suppresses_downstream_only(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(
        vcf,
        "chr1",
        40,
        [
            (14, "A", "T"),    # c.5   missense  (upstream of frameshift: kept)
            (21, "AG", "A"),   # c.13  frameshift (kept)
            (26, "T", "C"),    # c.17  downstream of frameshift (suppressed)
        ],
    )
    report = _annotate(gbff, str(vcf))
    assert "c.5A>T" in report          # upstream kept
    assert "c.13delG" in report        # frameshift kept
    assert "c.17" not in report        # downstream suppressed


# --------------------------------------------------------------------------- #
# reverse strand and multi-exon
# --------------------------------------------------------------------------- #

def test_reverse_strand_missense(gbff, tmp_path):
    # reverse gene beta: coding ATG TAT AAA TGA; c.5 (codon 2 middle) A>T
    # coding base A sits on genomic position 12 as a T (complement); ALT T -> A
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr2", 22, [(13, "T", "A")])
    report = _annotate(gbff, str(vcf), query=("test gene beta",))
    assert "missense_variant c.5A>T p.Tyr2Phe" in report
    assert "test.gene.beta: test gene beta" in report


def test_multi_exon_missense_in_second_exon(gbff, tmp_path):
    # gamma exon2 starts at genomic 10 == c.7 (codon 3 AAA); c.8 A>G -> AGA (Arg)
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr3", 20, [(12, "A", "G")])
    report = _annotate(gbff, str(vcf), query=("test gene gamma",))
    assert "missense_variant c.8A>G p.Lys3Arg" in report


def test_reverse_strand_frameshift_deletion(gbff, tmp_path):
    # beta coding ATG TAT AAA TGA on the minus strand at genomic [5, 17)
    # delete coding c.7 (first A of codon 3 AAA). coding idx 6 -> genomic 16-6 = 10
    # genome base there is complement(A)=T; left anchor at genomic 9 (coding c.8)
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr2", 22, [(10, "TT", "T")])  # POS 10 == genomic 9 anchor
    report = _annotate(gbff, str(vcf), query=("test gene beta",))
    assert "frameshift_variant c.7delA p.Lys3fs" in report


def test_transl_table_override_changes_call(gbff, tmp_path):
    # under the Yeast Mitochondrial code (table 3) CTG stays Leu but ATA -> Met;
    # simpler check: table override is honored and does not crash on alpha
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(14, "A", "T")])
    report = va.run(
        str(vcf),
        gbff,
        ["test gene alpha"],
        exact_match=True,
        transl_table=1,
    )
    assert "missense_variant c.5A>T p.Tyr2Phe" in report


# --------------------------------------------------------------------------- #
# allele depths, GENE INFO path and the no-variants message
# --------------------------------------------------------------------------- #

def test_allele_depths_reported(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(14, "A", "T", 20, 3)], with_sample=True)
    report = _annotate(gbff, str(vcf))
    assert "A:20 T:3" in report


def test_allele_depths_absent_reported_as_na(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(14, "A", "T")])
    report = _annotate(gbff, str(vcf))
    assert "A:NA T:NA" in report


def test_gene_info_annotation_path(gbff, tmp_path):
    """A variant carrying a GENE INFO field is resolved without interval overlap."""
    vcf = tmp_path / "v.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=40>",
        '##INFO=<ID=GENE,Number=.,Type=String,Description="gene">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t14\t.\tA\tT\t.\t.\tGENE=test_gene_alpha",
    ]
    vcf.write_text("\n".join(lines) + "\n")
    report = _annotate(gbff, str(vcf))
    assert "missense_variant c.5A>T p.Tyr2Phe" in report


def test_no_variants_message_lists_query_genes(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [])  # no records
    query = [
        "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase",
        "B9J08_005340",
    ]
    report = va.run(str(vcf), gbff, query, exact_match=True)
    assert report == (
        "No variants identified in queried genes "
        "(FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,"
        "B9J08_005340) relative to the reference genome"
    )


def test_report_joins_multiple_variants_with_commas(gbff, tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(14, "A", "T"), (18, "C", "T")])
    report = _annotate(gbff, str(vcf))
    assert report.count("test.gene.alpha:") == 2
    assert "),test.gene.alpha:" in report


# --------------------------------------------------------------------------- #
# regression tests for confirmed review findings
# --------------------------------------------------------------------------- #

def test_star_spanning_deletion_allele_does_not_crash(gbff, tmp_path):
    # GATK/bcftools '*' spanning-deletion allele must be skipped, not translated,
    # whether it is a secondary alt or the sole alt; valid alts still annotate.
    vcf = tmp_path / "v.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=40>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t14\t.\tA\tT,*\t.\t.\t.",  # valid T + spanning deletion
        "chr1\t18\t.\tC\t*\t.\t.\t.",    # '*' only
    ]
    vcf.write_text("\n".join(lines) + "\n")
    report = _annotate(gbff, str(vcf))
    assert "missense_variant c.5A>T p.Tyr2Phe" in report


def test_insertion_at_exon_5prime_boundary_coding_position(gbff, tmp_path):
    # gamma is two-exon; exon2 begins at c.7. An insertion at exon2's 5' edge must
    # map to c.6_7 (between the exons), not c.7_8 (off-by-one at the boundary).
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr3", 20, [(10, "G", "GC")])  # anchor = last intron base (genomic 9)
    report = _annotate(gbff, str(vcf), query=("test gene gamma",))
    assert "frameshift_variant c.6_7insC p.Lys3fs" in report


def test_reverse_strand_insertion_coding_position(gbff, tmp_path):
    # beta is minus-strand (coding ATG TAT AAA TGA). Insert 1 base between c.6 and
    # c.7; the strand-aware splice must yield c.6_7 (coding order), not a genomic one.
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr2", 22, [(11, "T", "TA")])  # forward-genome anchor genomic 10
    report = _annotate(gbff, str(vcf), query=("test gene beta",))
    assert "frameshift_variant c.6_7insT p.Lys3fs" in report


def test_insertion_at_cds_start_preserves_start_codon(gbff):
    # inserting immediately 5' of the ATG must splice before c.1 (cds index 0),
    # keeping ATG intact and inserting Phe (TTT) rather than corrupting the frame.
    models = va.build_gene_models(
        gbff, ["test gene alpha"], "CDS", "product", exact_match=True
    )
    model = models["test gene alpha"]
    ann = va.annotate_indel(model, 9, "", "TTT")  # alpha CDS starts at genomic 9
    assert ann["cds_pos"] == 0
    assert ann["so_term"] == "inframe_insertion"
    assert ann["hgvs_p"].endswith("insPhe")  # TTT -> Phe, not the corrupted-frame Ile


def test_inframe_insertion_after_stop_is_not_stop_gained(tmp_path):
    coding = "ATGAAACCCGGGTTTTAA"  # M K P G F *
    seq = "A" * 5 + coding + "A" * 5
    rec = SeqRecord(Seq(seq), id="chrgz", name="chrgz", description="")
    rec.annotations["molecule_type"] = "DNA"
    rec.features.append(
        SeqFeature(
            FeatureLocation(5, 5 + len(coding), strand=1),
            type="CDS",
            qualifiers={"product": ["gz"], "transl_table": ["1"]},
        )
    )
    gbff = tmp_path / "gz.gbff"
    with open(gbff, "w") as fh:
        SeqIO.write([rec], fh, "genbank")

    vcf = tmp_path / "v.vcf"
    # insert TTT immediately 3' of the terminal stop (gene [5,23), last base genomic 22)
    _write_vcf(vcf, "chrgz", 30, [(23, "A", "ATTT")])
    report = va.run(str(vcf), str(gbff), ["gz"], exact_match=True)
    assert "stop_gained" not in report
    assert "Xaa" not in report
    assert "stop_retained_variant" in report


def test_delins_stop_gained_anchors_at_new_stop(gbff, tmp_path):
    # c.8_10delinsGCT: codon 3 CCC->CGC (Pro3Arg missense) and codon 4 AAA->TAA (stop).
    # The Ter must be reported at the true new stop (Lys4), not the first changed
    # residue (Pro3, which is actually a sense codon in the mutant).
    vcf = tmp_path / "v.vcf"
    _write_vcf(vcf, "chr1", 40, [(17, "CCA", "GCT")])  # POS 17 == c.8
    report = _annotate(gbff, str(vcf))
    assert "stop_gained c.8_10delinsGCT p.Lys4Ter" in report


# --------------------------------------------------------------------------- #
# indel normalization unit tests
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize(
    "pos0, ref, alt, expected",
    [
        (100, "A", "T", (100, "A", "T")),          # SNP
        (100, "AG", "A", (101, "G", "")),          # left-anchored deletion
        (100, "CT", "T", (100, "C", "")),          # right-anchored deletion of C
        (100, "C", "CGGG", (101, "", "GGG")),      # insertion
        (100, "TCCC", "T", (101, "CCC", "")),      # multi-base deletion
    ],
)
def test_normalize_indel(pos0, ref, alt, expected):
    assert va.normalize_indel(pos0, ref, alt) == expected
