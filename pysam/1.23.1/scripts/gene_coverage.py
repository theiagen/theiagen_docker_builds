import sys
import json
import pysam
import logging
import argparse
from Bio import SeqIO
from collections import defaultdict


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def input_error_handling(args: argparse.ArgumentParser) -> None:
    """Handle incompatible input arguments"""
    if not args.bedfile and not args.reference_gbff:
        raise FileNotFoundError("'reference_gbff' or 'bedfile' is required")
    elif not args.query_genes and not args.bedfile:
        raise ValueError("'query_genes' or 'bedfile' required")


def exact_check(query_set: set, id: str) -> bool:
    """Return True or False for an exact match"""
    return id in query_set


def substring_check(query_set: set, id: str) -> bool:
    """Return True or False for a substring match"""
    return any(query in id for query in query_set)


def extract_queries_from_bed(bedfile):
    """Extract query regions from BED"""
    with open(bedfile, "r") as raw:
        return set(x.split()[3] for x in raw)


def write_json(filename: str, data: dict) -> None:
    """Write a JSON file compatible with WDL"""
    with open(filename, "w") as f:
        if data:
            json.dump(data, f, indent=4)
        else:
            # spoof Cromwell (Terra WDL)
            f.write('{"": 0}')


def parse_gbff(
    reference_gbff: str,
    query_set: set,
    feature_type: str,
    feature_qualifier: str,
    id_check: object,
    contig2query2coords: dict,
) -> dict:
    """Parse a GBFF to obtain query coordinates"""
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            record_id = record.name
            for feature in record.features:
                # is this the feature we want to scan?
                if feature.type.lower() == feature_type.lower():
                    # is there the qualifier that we want?
                    qualifier_ids = feature.qualifiers.get(feature_qualifier.strip())
                    if qualifier_ids:
                        qualifier_id = qualifier_ids[0]
                        # is this a qualifying feature?
                        if id_check(query_set, qualifier_id):
                            if qualifier_id in contig2query2coords[record_id]:
                                logger.warning(
                                    f"{qualifier_id} recovered multiple times in {record_id}"
                                )
                            # GenBanks are 1-based coordinates, though BioPython adjusts natively
                            loc_coords = [
                                [int(x.start), int(x.end)]
                                for x in feature.location.parts
                            ]
                            contig2query2coords[record_id][qualifier_id].extend(
                                loc_coords
                            )
    return contig2query2coords


def parse_bed(
    bedfile: str, query_set: set, id_check: object, contig2query2coords: dict
) -> dict:
    """Parse a BED file to obtain query coordinates"""
    with open(bedfile, "r") as handle:
        for line in handle:
            if not line.startswith("#"):
                data = line.split()
                id = data[3]
                # is this an entry we want?
                if id_check(query_set, id):
                    # BED files are 0-based coordinates
                    contig2query2coords[data[0]][id].append(
                        (int(data[1]), int(data[2]))
                    )
    return contig2query2coords


def import_bam(
    bamfile: str, contig2query2coords: dict, ambiguous_contig: bool
) -> tuple:
    imported_bam = pysam.AlignmentFile(bamfile)
    # generate an index if it does not exist
    if not imported_bam.has_index():
        logger.debug("Generating BAM index")
        pysam.index(bamfile)
        imported_bam = pysam.AlignmentFile(bamfile)

    # apply coordinates to first selected contig
    if ambiguous_contig:
        contig_names = imported_bam.references
        # can't apply ambiguous contig approach if there are multiple contigs
        if len(contig_names) > 1:
            raise ValueError(
                "can't use ambiguous_contig coordinates when there are multiple contigs in the reference"
            )
        contig = contig_names[0]
        # rename contig2query2coords to reflect first contig
        contig2query2coords = {contig: v for k, v in contig2query2coords.items()}
    return imported_bam, contig2query2coords


def quantify_gene_coverage(
    imported_bam: pysam.AlignmentFile,
    contig2query2coords: dict,
    min_depth: int = 1,
    min_quality: int = 1,
) -> tuple:
    """Quantify gene breadth and depth off coverage"""
    depth_dict = {}
    coverage_dict = {}
    reference_names = set(imported_bam.references)

    for contig, query2coords in contig2query2coords.items():
        if contig not in reference_names:
            raise ValueError(f"Contig '{contig}' not found in BAM references")
        contig_len = imported_bam.get_reference_length(contig)
        for query, loc_parts in query2coords.items():
            if query in depth_dict:
                logger.warning(
                    f"{query} is present on multiple contigs and will be overwritten"
                )
            # check coverage data across range
            depths = []
            coverages = []
            for coords in loc_parts:
                start, end = int(coords[0]), int(coords[1])
                if end <= start:
                    raise ValueError(
                        f"Invalid region for query '{query}' on contig '{contig}': start ({start}) must be < end ({end})"
                    )
                if start < 0:
                    raise ValueError(
                        f"Invalid region for query '{query}' on contig '{contig}': start ({start}) must be >= 0"
                    )
                if end > contig_len:
                    raise ValueError(
                        f"Invalid region for query '{query}' on contig '{contig}': end ({end}) exceeds contig length ({contig_len})"
                    )
                coverage_data = imported_bam.count_coverage(
                    contig, start, end, quality_threshold=min_quality
                )
                for i, _ in enumerate(range(start, end)):
                    # calculate total depth across bases
                    total_depth = (
                        coverage_data[0][i]
                        + coverage_data[1][i]
                        + coverage_data[2][i]
                        + coverage_data[3][i]
                    )
                    # base is considered covered if beyond minimum depth
                    coverages.append(total_depth >= min_depth)
                    depths.append(total_depth)
            if not depths:
                raise ValueError(
                    f"No positions evaluated for query '{query}' on contig '{contig}'"
                )
            depth_dict[query] = sum(depths) / len(depths)
            # breadth is percent of covered bases exceeding min_depth
            coverage_dict[query] = 100 * (sum(coverages) / len(coverages))

    return dict(sorted(depth_dict.items())), dict(sorted(coverage_dict.items()))


def make_tsv(depth_dict: dict, coverage_dict: dict, ambiguous_contig: bool) -> str:
    """Make a readable TSV to convey depth and coverage"""
    if ambiguous_contig:
        name = "query (WARNING: results may be inaccurate if sample is not mapped to reference used to generate BED file coordinates)"
    else:
        name = "query"
    tsv_str = f"#{name}\taverage_depth\tpercent_coverage\n"
    for query, depth in depth_dict.items():
        tsv_str += f"{query}\t{depth}\t{coverage_dict[query]}\n"
    return tsv_str.strip()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="task_gene_coverage.wdl dependency script (github.com/theiagen/public_health_bioinformatics)"
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bedfile")
    parser.add_argument("--reference_gbff")
    parser.add_argument("--query_genes", nargs="+")
    parser.add_argument("--feature_type", default="CDS")
    parser.add_argument("--feature_qualifier", default="product")
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument("--ambiguous_contig", action="store_true")
    parser.add_argument("--min_depth", type=int, default=1)
    parser.add_argument("--min_quality", type=int, default=0)
    args = parser.parse_args()

    # error parsing
    input_error_handling(args)

    # set comparison check function
    if args.exact_match:
        id_check = exact_check
    else:
        id_check = substring_check

    # import queries
    if args.query_genes:
        query_set = set()
        for queries in args.query_genes:
            query_set = query_set.union(q.strip() for q in queries.split(","))
    else:
        query_set = extract_queries_from_bed(args.bedfile)

    # {<CONTIG>: <QUERY>: [(LOC_START_1, LOC_END_1,), (LOC_START_n, LOC_END_n),]}
    contig2query2coords = defaultdict(lambda: defaultdict(list))
    if args.reference_gbff:
        contig2query2coords = parse_gbff(
            args.reference_gbff,
            query_set,
            args.feature_type,
            args.feature_qualifier,
            id_check,
            contig2query2coords,
        )
    if args.bedfile:
        contig2query2coords = parse_bed(
            args.bedfile, query_set, id_check, contig2query2coords
        )

    # import BAM and modify contig coordinates if needed
    imported_bam, contig2query2coords = import_bam(
        args.bam, contig2query2coords, args.ambiguous_contig
    )

    # quantify statistics and write
    depth_dict, coverage_dict = quantify_gene_coverage(
        imported_bam, contig2query2coords, args.min_depth, args.min_quality
    )
    write_json("DEPTH_DICT.json", depth_dict)
    write_json("COVERAGE_DICT.json", coverage_dict)

    # add missing entries to TSV report
    missing_genes = query_set.difference(set(coverage_dict.keys()))
    for gene in missing_genes:
        # depth may be reported for those that have no breadth
        if gene not in depth_dict:
            depth_dict[gene] = 0
        coverage_dict[gene] = 0

    tsv_str = make_tsv(
        depth_dict, coverage_dict, args.ambiguous_contig and args.bedfile
    )
    with open("COVERAGE_STATS.tsv", "w") as out:
        out.write(tsv_str)

    sys.exit(0)
