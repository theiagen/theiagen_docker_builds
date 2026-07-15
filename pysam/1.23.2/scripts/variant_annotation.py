"""Annotate the protein-level consequences of variants that overlap query genes.

This module is the companion to ``gene_coverage.py``.  ``gene_coverage.py``
extracts the variants that overlap a set of query genes into ``GENE_VARIANTS.vcf``
(annotating each record with a ``GENE`` INFO field).  This module consumes that
VCF, together with the reference GenBank (GBFF) that supplies the coding
sequence, strand, product name and translation table for each gene, and reports
the effect of every variant on the translated protein.

For each variant it determines whether a substitution is

* ``synonymous_variant`` (silent, same amino acid),
* ``missense_variant``   (different amino acid),
* ``stop_gained``        (nonsense, a new stop codon),
* ``stop_lost``/``start_lost`` (edge substitutions),

and whether an indel is an ``inframe_deletion``, ``inframe_insertion`` or a
``frameshift_variant``.  A frameshift invalidates every downstream codon, so all
variants that fall 3' of a frameshift within the same gene are disregarded.

The report is emitted as a single string, one entry per surviving variant::

    <gene_id>: <product> (<so_term> <hgvs_c> <hgvs_p>; <REF>:<depth> <ALT>:<depth>)

with entries joined by commas.  When no variants are recovered a single sentence
naming the queried genes is emitted instead.
"""

import sys
import logging
import argparse
from collections import defaultdict

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# single-letter -> three-letter amino acid codes, with HGVS extensions
_AA_1TO3 = {k.upper(): v for k, v in IUPACData.protein_letters_1to3.items()}
_AA_1TO3["*"] = "Ter"
_AA_1TO3["X"] = "Xaa"

_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

_VALID_BASES = set("ACGTNacgtn")


def is_nucleotide_allele(allele) -> bool:
    """True only for plain nucleotide alleles.

    Filters symbolic alleles ('<DEL>', '<*>'), the spanning-deletion allele
    ('*') emitted by GATK/bcftools, and empty/missing alleles, all of which
    must not be translated as substitutions."""
    return bool(allele) and all(base in _VALID_BASES for base in allele)


def aa3(aa: str) -> str:
    """Return the three-letter code for a one-letter amino acid symbol"""
    return _AA_1TO3.get(aa.upper(), "Xaa")


def complement(seq: str) -> str:
    """Complement (not reverse) a nucleotide string"""
    return seq.translate(_COMPLEMENT)


def normalize_name(name: str) -> str:
    """Collapse whitespace and VCF-reserved characters to '.' for a stable id.

    e.g. 'lanosterol 14-alpha demethylase' -> 'lanosterol.14-alpha.demethylase'"""
    out = name
    for char in (" ", "\t", ";", "=", ","):
        out = out.replace(char, ".")
    # collapse runs of dots that can arise from adjacent replaced characters
    while ".." in out:
        out = out.replace("..", ".")
    return out.strip(".")


def sanitize_info_value(value: str) -> str:
    """Sanitize a string for use in a VCF INFO field (mirrors gene_coverage.py)"""
    for char in (" ", "\t", ";", "=", ","):
        value = value.replace(char, "_")
    return value


class GeneModel:
    """A coding sequence model for a single query gene.

    Coordinates are 0-based, half-open and refer to the reference contig.
    ``genomic_positions`` lists every coding base position in translation
    (5'->3') order; ``ref_coding`` is the reference coding sequence in that same
    order (reverse-complemented for minus-strand genes); ``pos2cds`` maps a
    genomic position to its index within ``ref_coding``."""

    def __init__(self, gene_id, product, contig, strand, transl_table):
        self.gene_id = gene_id
        self.product = product
        self.contig = contig
        self.strand = strand
        self.transl_table = transl_table
        self.genomic_positions = []
        self.pos2cds = {}
        self.ref_coding = ""
        self.ref_protein = ""
        # genomic span used for interval-overlap gene assignment
        self.genomic_start = None
        self.genomic_end = None

    def finalize(self, contig_seq: str) -> None:
        """Build the coding sequence, position index and reference protein"""
        bases = [contig_seq[p].upper() for p in self.genomic_positions]
        coding = "".join(bases)
        if self.strand == -1:
            # positions are already reversed; complement each base to get revcomp
            coding = complement(coding)
        self.ref_coding = coding
        self.pos2cds = {g: i for i, g in enumerate(self.genomic_positions)}
        self.ref_protein = translate(coding, self.transl_table)
        if self.genomic_positions:
            self.genomic_start = min(self.genomic_positions)
            self.genomic_end = max(self.genomic_positions) + 1

    def codon(self, codon_number: int) -> str:
        """Return the reference codon (1-based codon number) or '' if incomplete"""
        start = (codon_number - 1) * 3
        codon = self.ref_coding[start : start + 3]
        return codon if len(codon) == 3 else ""

    def aa_at(self, codon_number: int) -> str:
        """Return the reference amino acid (one letter) at a 1-based codon"""
        codon = self.codon(codon_number)
        return translate(codon, self.transl_table) if codon else "X"


def translate(seq: str, table) -> str:
    """Translate a nucleotide string, truncating to a whole number of codons"""
    trimmed = seq[: len(seq) - (len(seq) % 3)]
    if not trimmed:
        return ""
    return str(Seq(trimmed).translate(table=table))


def _ordered_genomic_positions(parts, strand):
    """Ordered 0-based genomic positions of a CDS in translation (5'->3') order"""
    ordered = sorted(((int(s), int(e)) for s, e in parts), key=lambda x: x[0])
    positions = []
    for start, end in ordered:
        positions.extend(range(start, end))
    if strand == -1:
        positions.reverse()
    return positions


def _match_query(query_list, identifiers, exact_match: bool):
    """Return the first query term (in order) matching any candidate identifier.

    Matching is normalization-aware, so dotted query names (e.g.
    'lanosterol.14-alpha.demethylase') match spaced GBFF products
    ('lanosterol 14-alpha demethylase') and vice versa."""
    for query in query_list:
        nq = normalize_name(query)
        for ident in identifiers:
            ni = normalize_name(ident)
            if exact_match:
                if query == ident or nq == ni:
                    return query
            elif query in ident or nq in ni:
                return query
    return None


def build_gene_models(
    reference_gbff: str,
    query_genes,
    feature_type: str,
    feature_qualifier: str,
    exact_match: bool = False,
    transl_table_override: int = None,
) -> dict:
    """Parse a GBFF into {lookup_key: GeneModel} for every query gene.

    A CDS is matched by the ``feature_qualifier`` value plus its ``gene``,
    ``locus_tag`` and ``protein_id`` qualifiers, so query sets may mix product
    names and locus tags.  A model is registered under many lookup keys (raw,
    sanitized and normalized forms of every identifier) so it can be recovered
    from whatever identifier the VCF ``GENE`` field carries."""
    query_list = list(query_genes)
    id_qualifiers = [feature_qualifier.strip(), "gene", "locus_tag", "protein_id"]
    models_by_key = {}
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_seq = str(record.seq)
            for feature in record.features:
                if feature.type.lower() != feature_type.lower():
                    continue
                identifiers = []
                for qualifier in id_qualifiers:
                    identifiers.extend(feature.qualifiers.get(qualifier, []))
                if not identifiers:
                    continue

                matched_query = _match_query(query_list, identifiers, exact_match)
                if matched_query is None:
                    continue

                strand = feature.location.strand
                if strand not in (1, -1):
                    logger.warning(
                        f"Skipping '{matched_query}' on {record.name}: "
                        f"unresolved strand ({strand})"
                    )
                    continue

                if transl_table_override is not None:
                    transl_table = transl_table_override
                else:
                    transl_table = int(
                        feature.qualifiers.get("transl_table", ["1"])[0]
                    )

                product_vals = feature.qualifiers.get(feature_qualifier.strip())
                product = product_vals[0] if product_vals else matched_query
                gene_id = normalize_name(matched_query)
                model = GeneModel(
                    gene_id=gene_id,
                    product=product,
                    contig=record.name,
                    strand=strand,
                    transl_table=transl_table,
                )
                parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
                model.genomic_positions = _ordered_genomic_positions(parts, strand)
                model.finalize(contig_seq)

                keys = set()
                for ident in identifiers + [matched_query, product, gene_id]:
                    keys.update(
                        (ident, sanitize_info_value(ident), normalize_name(ident))
                    )
                for key in keys:
                    if not key:
                        continue
                    if key in models_by_key and models_by_key[key] is not model:
                        logger.warning(
                            f"'{key}' recovered multiple times; keeping first"
                        )
                    else:
                        models_by_key[key] = model
    return models_by_key


def normalize_indel(pos0: int, ref: str, alt: str) -> tuple:
    """Trim shared prefix/suffix from a REF/ALT pair.

    Returns (changed_pos0, ref_segment, alt_segment) where ``changed_pos0`` is
    the 0-based genomic position of the first base of ``ref_segment``.  For a
    pure insertion ``ref_segment`` is empty and the insertion falls between
    ``changed_pos0 - 1`` and ``changed_pos0``; for a pure deletion
    ``alt_segment`` is empty."""
    ref = ref.upper()
    alt = alt.upper()
    # trim shared prefix
    prefix = 0
    while prefix < len(ref) and prefix < len(alt) and ref[prefix] == alt[prefix]:
        prefix += 1
    # trim shared suffix (not past the trimmed prefix)
    suffix = 0
    while (
        suffix < len(ref) - prefix
        and suffix < len(alt) - prefix
        and ref[len(ref) - 1 - suffix] == alt[len(alt) - 1 - suffix]
    ):
        suffix += 1
    ref_seg = ref[prefix : len(ref) - suffix]
    alt_seg = alt[prefix : len(alt) - suffix]
    return pos0 + prefix, ref_seg, alt_seg


def allele_depths(record, alt_index: int):
    """Return (ref_depth, alt_depth) from the first sample's AD field, or (None, None)"""
    if not len(record.samples):
        return None, None
    sample = next(iter(record.samples.values()))
    try:
        ad = sample["AD"]
    except (KeyError, ValueError, TypeError):
        ad = None
    if ad is None:
        return None, None
    ref_depth = ad[0] if len(ad) > 0 else None
    alt_depth = ad[alt_index + 1] if len(ad) > alt_index + 1 else None
    return ref_depth, alt_depth


def _depth_str(depth) -> str:
    return "NA" if depth is None else str(depth)


def annotate_snp(model: GeneModel, pos0: int, ref: str, alt: str) -> dict:
    """Annotate a single-nucleotide substitution against a gene's coding sequence"""
    cds_idx = model.pos2cds.get(pos0)
    if cds_idx is None:
        return None  # not within a coding base (e.g. intronic anchor)

    coding_ref = ref if model.strand == 1 else complement(ref)
    coding_alt = alt if model.strand == 1 else complement(alt)

    expected = model.ref_coding[cds_idx]
    if expected != coding_ref:
        logger.warning(
            f"{model.gene_id}: reference base '{expected}' at c.{cds_idx + 1} "
            f"disagrees with VCF REF '{coding_ref}' (coding strand)"
        )

    codon_number = cds_idx // 3 + 1
    pos_in_codon = cds_idx % 3
    ref_codon = model.codon(codon_number)
    hgvs_c = f"c.{cds_idx + 1}{coding_ref}>{coding_alt}"

    if not ref_codon:
        # incomplete terminal codon; report nucleotide change only
        return {
            "so_term": "coding_sequence_variant",
            "hgvs_c": hgvs_c,
            "hgvs_p": "p.?",
            "cds_pos": cds_idx,
            "is_frameshift": False,
        }

    mut_codon = ref_codon[:pos_in_codon] + coding_alt + ref_codon[pos_in_codon + 1 :]
    ref_aa = translate(ref_codon, model.transl_table)
    alt_aa = translate(mut_codon, model.transl_table)

    if ref_aa == alt_aa:
        so_term = "synonymous_variant"
        hgvs_p = f"p.{aa3(ref_aa)}{codon_number}="
    elif alt_aa == "*":
        so_term = "stop_gained"
        hgvs_p = f"p.{aa3(ref_aa)}{codon_number}Ter"
    elif ref_aa == "*":
        so_term = "stop_lost"
        hgvs_p = f"p.Ter{codon_number}{aa3(alt_aa)}"
    elif codon_number == 1 and ref_aa == "M":
        so_term = "start_lost"
        hgvs_p = f"p.{aa3(ref_aa)}1{aa3(alt_aa)}"
    else:
        so_term = "missense_variant"
        hgvs_p = f"p.{aa3(ref_aa)}{codon_number}{aa3(alt_aa)}"

    return {
        "so_term": so_term,
        "hgvs_c": hgvs_c,
        "hgvs_p": hgvs_p,
        "cds_pos": cds_idx,
        "is_frameshift": False,
    }


def annotate_indel(
    model: GeneModel, changed_pos0: int, ref_seg: str, alt_seg: str
) -> dict:
    """Annotate an insertion, deletion or delins against a gene's coding sequence"""
    # coding indices of the deleted reference bases (those that are coding)
    deleted_idx = sorted(
        model.pos2cds[changed_pos0 + i]
        for i in range(len(ref_seg))
        if (changed_pos0 + i) in model.pos2cds
    )
    coding_alt = alt_seg if model.strand == 1 else complement(alt_seg)
    if model.strand == -1:
        coding_alt = coding_alt[::-1]  # reverse to restore 5'->3' coding order

    deleted_cds = len(deleted_idx)
    inserted_cds = len(coding_alt)
    frame_delta = (inserted_cds - deleted_cds) % 3

    # locate the variant within the coding sequence for splicing / reporting
    if deleted_idx:
        cds_start = deleted_idx[0]
        cds_end = deleted_idx[-1] + 1
    else:
        # pure insertion: require a coding base immediately on one side, then
        # splice at the count of coding bases lying 5' of the insertion point.
        # This is strand-aware and correct at exon/CDS boundaries where a naive
        # flank index would be off by one (5' base of an internal exon / CDS start)
        if (changed_pos0 - 1) not in model.pos2cds and changed_pos0 not in model.pos2cds:
            return None  # insertion is not adjacent to any coding base
        if model.strand == 1:
            cds_start = sum(1 for g in model.genomic_positions if g < changed_pos0)
        else:
            cds_start = sum(1 for g in model.genomic_positions if g >= changed_pos0)
        cds_end = cds_start

    first_codon = cds_start // 3 + 1

    if frame_delta != 0:
        return {
            "so_term": "frameshift_variant",
            "hgvs_c": _hgvs_c_indel(cds_start, cds_end, coding_alt, ref_seg, model),
            "hgvs_p": f"p.{aa3(model.aa_at(first_codon))}{first_codon}fs",
            "cds_pos": cds_start,
            "is_frameshift": True,
        }

    # in-frame: build the mutant coding sequence and compare proteins
    mut_coding = model.ref_coding[:cds_start] + coding_alt + model.ref_coding[cds_end:]
    mut_protein = translate(mut_coding, model.transl_table)
    ref_protein = model.ref_protein

    hgvs_c = _hgvs_c_indel(cds_start, cds_end, coding_alt, ref_seg, model)
    so_term, hgvs_p = _protein_consequence(
        ref_protein, mut_protein, first_codon, deleted_cds, inserted_cds
    )
    return {
        "so_term": so_term,
        "hgvs_c": hgvs_c,
        "hgvs_p": hgvs_p,
        "cds_pos": cds_start,
        "is_frameshift": False,
    }


def _hgvs_c_indel(cds_start, cds_end, coding_alt, ref_seg, model) -> str:
    """Build the HGVS c. string for an indel/delins in coding coordinates"""
    deleted_len = cds_end - cds_start
    if not ref_seg or deleted_len == 0:
        # pure insertion between cds_start-1 and cds_start (1-based cds_start, +1)
        return f"c.{cds_start}_{cds_start + 1}ins{coding_alt}"
    deleted_bases = model.ref_coding[cds_start:cds_end]
    if not coding_alt:
        if deleted_len == 1:
            return f"c.{cds_start + 1}del{deleted_bases}"
        return f"c.{cds_start + 1}_{cds_end}del"
    # delins
    if deleted_len == 1:
        return f"c.{cds_start + 1}delins{coding_alt}"
    return f"c.{cds_start + 1}_{cds_end}delins{coding_alt}"


def _protein_consequence(
    ref_protein, mut_protein, first_codon, deleted_cds, inserted_cds
) -> tuple:
    """Classify an in-frame protein change and build the HGVS p. string"""
    net = inserted_cds - deleted_cds
    net_aa = net // 3  # frame_delta == 0, so net is a whole number of codons

    # translated products up to (and excluding) the first stop codon
    ref_prod = ref_protein.split("*", 1)[0]
    mut_prod = mut_protein.split("*", 1)[0]
    first_res = first_codon - 1  # 0-based first affected residue

    # a variant lying entirely 3' of the reference stop codon cannot alter the
    # protein product; only report a change if it actually reaches the terminator
    if first_res >= len(ref_prod):
        if mut_prod == ref_prod:
            return "stop_retained_variant", "p.(=)"
        return "stop_lost", f"p.Ter{len(ref_prod) + 1}ext"

    # a genuine stop_gained truncates the product earlier than the indel's own
    # length change accounts for; a plain in-frame deletion only shortens it by
    # net_aa residues and must not be mistaken for a new stop. Anchor the Ter at
    # the mutant's first stop, not at the first *changed* residue (which may be a
    # sense codon upstream of the new stop in a multi-codon delins)
    if len(mut_prod) < len(ref_prod) + net_aa:
        stop_pos = mut_protein.find("*")
        codon = stop_pos + 1
        ref_res = _res(ref_protein, stop_pos - net_aa)
        return "stop_gained", f"p.{aa3(ref_res)}{codon}Ter"

    idx = _first_diff(ref_protein, mut_protein)

    if net < 0:  # net deletion of amino acids
        n_del = (-net) // 3
        first = idx
        last = idx + n_del - 1
        if n_del <= 1:
            hgvs_p = f"p.{aa3(_res(ref_protein, first))}{first + 1}del"
        else:
            hgvs_p = (
                f"p.{aa3(_res(ref_protein, first))}{first + 1}_"
                f"{aa3(_res(ref_protein, last))}{last + 1}del"
            )
        return "inframe_deletion", hgvs_p

    if net > 0:  # net insertion of amino acids
        n_ins = net // 3
        before = idx - 1
        after = idx
        inserted = "".join(aa3(_res(mut_protein, idx + k)) for k in range(n_ins))
        hgvs_p = (
            f"p.{aa3(_res(ref_protein, before))}{before + 1}_"
            f"{aa3(_res(ref_protein, after))}{after + 1}ins{inserted}"
        )
        return "inframe_insertion", hgvs_p

    # equal length: one or more codons substituted (MNV / in-frame delins)
    if idx >= len(ref_protein) or idx >= len(mut_protein):
        return "synonymous_variant", "p.="
    if _res(ref_protein, idx) == _res(mut_protein, idx):
        return "synonymous_variant", "p.="
    return (
        "missense_variant",
        f"p.{aa3(_res(ref_protein, idx))}{idx + 1}{aa3(_res(mut_protein, idx))}",
    )


def _first_diff(a: str, b: str) -> int:
    """Index of the first position at which two strings differ"""
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return i
    return min(len(a), len(b))


def _res(protein: str, idx: int) -> str:
    """Residue at ``idx`` (one letter), or 'X' if out of range"""
    return protein[idx] if 0 <= idx < len(protein) else "X"


def genes_for_record(record, models_by_key: dict, interval_index: dict) -> list:
    """Resolve the GeneModel(s) a record belongs to.

    Prefers the VCF ``GENE`` INFO annotation; falls back to interval overlap so
    the module also works on a raw (un-extracted) VCF."""
    models = []
    seen = set()
    gene_info = record.info.get("GENE") if "GENE" in record.info else None
    if gene_info:
        names = gene_info if isinstance(gene_info, (list, tuple)) else [gene_info]
        for name in names:
            model = models_by_key.get(name)
            if model is not None and id(model) not in seen:
                seen.add(id(model))
                models.append(model)
    if not models:
        for start, end, model in interval_index.get(record.contig, []):
            if record.start < end and record.stop > start and id(model) not in seen:
                seen.add(id(model))
                models.append(model)
    return models


def annotate_vcf(vcffile: str, models_by_key: dict) -> list:
    """Annotate every query-gene variant in a VCF.

    Returns a list of annotation dicts in VCF read order, after dropping every
    variant that lies downstream (3') of a frameshift within the same gene."""
    interval_index = defaultdict(list)
    for model in set(models_by_key.values()):
        if model.genomic_start is not None:
            interval_index[model.contig].append(
                (model.genomic_start, model.genomic_end, model)
            )

    annotations = []
    read_order = 0
    with pysam.VariantFile(vcffile) as vcf_in:
        for record in vcf_in:
            models = genes_for_record(record, models_by_key, interval_index)
            if not models or not record.alts or not is_nucleotide_allele(record.ref):
                continue
            for alt_index, alt in enumerate(record.alts):
                # skip symbolic ('<...>'), spanning-deletion ('*') and NON-REF alleles
                if not is_nucleotide_allele(alt):
                    continue
                ref_depth, alt_depth = allele_depths(record, alt_index)
                changed_pos0, ref_seg, alt_seg = normalize_indel(
                    record.start, record.ref, alt
                )
                for model in models:
                    try:
                        if len(ref_seg) == 1 and len(alt_seg) == 1:
                            ann = annotate_snp(model, changed_pos0, ref_seg, alt_seg)
                        elif not ref_seg and not alt_seg:
                            ann = None
                        else:
                            ann = annotate_indel(model, changed_pos0, ref_seg, alt_seg)
                    except Exception as exc:  # never let one record abort the report
                        logger.warning(
                            f"Skipping {record.contig}:{record.pos} "
                            f"{record.ref}>{alt} for {model.gene_id}: {exc}"
                        )
                        ann = None
                    if ann is None:
                        continue
                    ann.update(
                        {
                            "gene_id": model.gene_id,
                            "product": model.product,
                            "ref_allele": record.ref,
                            "alt_allele": alt,
                            "ref_depth": ref_depth,
                            "alt_depth": alt_depth,
                            "read_order": read_order,
                        }
                    )
                    annotations.append(ann)
            read_order += 1

    return _apply_frameshift_suppression(annotations)


def _apply_frameshift_suppression(annotations: list) -> list:
    """Drop variants that fall 3' of the first frameshift within each gene"""
    first_fs_pos = {}
    for ann in annotations:
        if ann["is_frameshift"]:
            gene = ann["gene_id"]
            pos = ann["cds_pos"]
            if gene not in first_fs_pos or pos < first_fs_pos[gene]:
                first_fs_pos[gene] = pos

    kept = []
    for ann in annotations:
        cutoff = first_fs_pos.get(ann["gene_id"])
        if cutoff is not None and ann["cds_pos"] > cutoff:
            logger.debug(
                f"Suppressing {ann['gene_id']} {ann['hgvs_c']} "
                f"(downstream of frameshift at c.{cutoff + 1})"
            )
            continue
        kept.append(ann)

    kept.sort(key=lambda a: a["read_order"])
    return kept


def format_report(annotations: list, ordered_query_genes) -> str:
    """Render the annotation report string (or the no-variants sentence)"""
    if not annotations:
        genes = ",".join(ordered_query_genes)
        return (
            f"No variants identified in queried genes ({genes}) "
            "relative to the reference genome"
        )

    entries = []
    for ann in annotations:
        depths = (
            f"{ann['ref_allele']}:{_depth_str(ann['ref_depth'])} "
            f"{ann['alt_allele']}:{_depth_str(ann['alt_depth'])}"
        )
        entries.append(
            f"{ann['gene_id']}: {ann['product']} "
            f"({ann['so_term']} {ann['hgvs_c']} {ann['hgvs_p']}; {depths})"
        )
    return ",".join(entries)


def ordered_query_genes(query_genes_arg) -> list:
    """Flatten the --query_genes argument into an ordered, de-duplicated list"""
    ordered = []
    seen = set()
    for chunk in query_genes_arg or []:
        for gene in chunk.split(","):
            gene = gene.strip()
            if gene and gene not in seen:
                seen.add(gene)
                ordered.append(gene)
    return ordered


def run(
    vcffile: str,
    reference_gbff: str,
    query_genes,
    feature_type: str = "CDS",
    feature_qualifier: str = "product",
    exact_match: bool = False,
    transl_table: int = None,
) -> str:
    """Annotate a VCF and return the report string"""
    query_arg = query_genes if isinstance(query_genes, (list, tuple)) else [query_genes]
    ordered = ordered_query_genes(query_arg)
    models_by_key = build_gene_models(
        reference_gbff,
        ordered,
        feature_type,
        feature_qualifier,
        exact_match=exact_match,
        transl_table_override=transl_table,
    )
    annotations = annotate_vcf(vcffile, models_by_key)
    return format_report(annotations, ordered)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate protein-level consequences of gene-overlapping variants "
        "(task_gene_coverage.wdl dependency; github.com/theiagen/public_health_bioinformatics)"
    )
    parser.add_argument("--vcf", required=True, help="VCF of gene-overlapping variants")
    parser.add_argument(
        "--reference_gbff",
        required=True,
        help="reference GenBank supplying coding sequence, strand and product",
    )
    parser.add_argument("--query_genes", nargs="+", required=True)
    parser.add_argument("--feature_type", default="CDS")
    parser.add_argument("--feature_qualifier", default="product")
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument(
        "--transl_table",
        type=int,
        default=None,
        help="override the genetic code (default: read /transl_table from each CDS)",
    )
    parser.add_argument("--output", default="VARIANT_ANNOTATIONS.txt")
    args = parser.parse_args()

    report = run(
        args.vcf,
        args.reference_gbff,
        args.query_genes,
        feature_type=args.feature_type,
        feature_qualifier=args.feature_qualifier,
        exact_match=args.exact_match,
        transl_table=args.transl_table,
    )
    with open(args.output, "w") as out:
        out.write(report + "\n")
    print(report)
    sys.exit(0)
