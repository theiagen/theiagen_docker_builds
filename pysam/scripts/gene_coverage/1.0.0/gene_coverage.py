import sys
import pysam
import logging
import argparse
from Bio import SeqIO
from collections import defaultdict


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def exact_check(query_set: set, id: str) -> bool:
    """Return True or False for an exact match"""
    return id in query_set


def substring_check(query_set: set, id: str) -> bool:
    """Return True or False for a substring match"""
    return any(id in query for query in query_set)


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
    id_check: function,
    contig2query2coords: dict,
) -> dict:
    """Parse a GBFF to obtain query coordinates"""
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                # is this the feature we want to scan?
                if feature.type.lower() == feature_type:
                    # is there the qualifier that we want?
                    qualifier_id = feature.qualifiers.get(feature_qualifier).strip()
                    if qualifier_id:
                        # is this a qualifying feature?
                        if id_check(query_set, qualifier_id):
                            if qualifier_id in contig2query2coords[record.name]:
                                print(
                                    f"WARNING: {qualifier_id} recovered multiple times in {record.name}"
                                )
                            # GenBanks are 1-based coordinates, which we need to adjust for PySam
                            contig2query2coords[record.name][qualifier_id] = (
                                int(feature.location.start) - 1,
                                int(feature.location.end) - 1,
                            )
    return contig2query2coords


def parse_bed(
    bedfile: str, query_set: set, id_check: function, contig2query2coords: dict
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
                    contig2query2coords[data[0]][id] = (int(data[1]), int(data[2]))
    return contig2query2coords


def quantify_gene_coverage(
    bamfile: str, contig2query2coords: dict, min_depth: int = 1
) -> tuple:
    """Quantify gene breadth and depth off coverage"""
    imported_bam = pysam.ALignmentFile(bamfile)

    depth_dict = {}
    coverage_dict = {}

    for contig in query2coords in contig2query2coords.items():
        for query, coords in query2coords.items():
            if query in depth_dict:
                logger.warning(
                    f"{query} is present on multiple contigs and will be overwritten"
                )
            # check coverage data across range
            coverage_data = imported_bam.count_coverage(contig, coords[0], coords[1])
            depths = []
            coverages = []
            for i, pos in enumerate(range(coords[0], coords[1])):
                # calculate total depth across bases
                total_depth = sum(
                    coverage_data[0][i]
                    + coverage_data[1][i]
                    + coverage_data[2][i]
                    + coverage_data[3][i]
                )
                # base is considered covered if beyond minimum depth
                coverages.append(total_depth >= min_depth)
                depths.append(total_depth)
            depth_dict[query] = sum(depths) / len(depths)
            # breadth is percent of covered bases exceeding min_depth
            coverage_dict[query] = 100 * (sum(coverages) / len(coverages))

    return depth_dict, coverage_dict


def make_csv(depth_dict: dict, coverage_dict: dict) -> str:
    """Make a readable CSV to convey depth and coverage"""
    csv_str = "#gene,average_depth,percent_coverage\n"
    for query, depth in depth_dict.items():
        csv_str += f"{query},{depth},{coverage_dict[query]}\n"
    return csv_str.strip()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="task_gene_coverage.wdl dependency script (github.com/theiagen/public_health_bioinformatics)"
    )
    parser.add_argument("--bamfile", required=True)
    parser.add_argument("--bedfile")
    parser.add_argument("--reference_gbff")
    parser.add_argument("--query_genes", required=True, nargs="+")
    parser.add_argument("--feature_type", default="mRNA")
    parser.add_argument("--feature_qualifier", default="product")
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument("--min_depth", type=int, default=1)
    args = parser.parse_args()

    if not args.bedfile and not args.reference_gbff:
        raise FileNotFoundError("'reference_gbff' or 'bedfile' is required")

    # set comparison check function
    if args.exact_match:
        id_check = exact_check
    else:
        id_check = substring_check

    # import queries
    query_set = set()
    for queries in args.query_genes:
        query_set = query_set.union(q.strip() for q in queries.split(","))

    contig2query2coords = defaultdict(dict)
    if args.reference_gbff:
        contig2query2coords = parse_gbff(
            args.reference_gbff,
            query_set,
            args.feature_type,
            args.feature_qualifier,
            id_check,
            contig2query2coords,
            id_check,
        )
    if args.bedfile:
        contig2query2coords = parse_bed(
            args.bedfile, query_set, id_check, contig2query2coords
        )

    depth_dict, coverage_dict = quantify_gene_coverage(
        args.bamfile, contig2query2coords, args.min_depth
    )

    csv_str = make_csv(depth_dict, coverage_dict)
    with open("COVERAGE_STATS.csv", "w") as out:
        out.write(csv_str)


    write_json("DEPTH_DICT.json", depth_dict)
    write_json("COVERAGE_DICT.json", coverage_dict)

    sys.exit(0)
