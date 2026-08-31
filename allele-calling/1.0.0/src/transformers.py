from collections import Counter, defaultdict
from functools import partial
from logging import Logger
from math import floor
from pathlib import Path
from tempfile import NamedTemporaryFile
from xml.etree.ElementTree import fromstring

import numpy as np
from ngs_pipeline_lib.base.inputs import FastaInput, FilePath
from ngs_pipeline_lib.tools.parallel import run_parallel
from pysam import AlignedSegment, AlignmentFile, PileupRead
from pysam import index as pysam_index
from pysam import sort as pysam_sort
from pysam import view as pysam_view

from src.models import AlleleCall, FlaggedAlleleCall
from src.sam_flags import CombinedFlag, Flag
from src.tools import LocusGroups, get_contigs_seq

BASIC_NTS = {"A", "T", "C", "G"}
STRANDS = {"Fwd", "Rev"}
MULTIPLE_CALLS_FLAGS = [Flag.DOUBLE_CALLS, Flag.MORE_CALLS]


class XMLdata:
    def __init__(
        self,
        xml: Path,
    ):
        self.calls: dict[str, list[AlleleCall]] = {}

        with open(xml, encoding="utf-8") as reader:
            self.parse_xml(reader.read())

    def parse_xml(
        self,
        xml_content: str,
    ):
        """
        Structure of the xml:
        <loci>
            <l>
                <id>SALM_25362</id>                                         Locus identifier
                <as>                                                        A tag with one or more <a>
                <a>                                                         Contains all the info for a single allele call
                        <id val="7" />                                      Allele identifier
                        <si type="valid" val="100" />                       Sequence identity, equivalent to similarity
                        <ev type="valid" val="0" />                         E-value (from blast)
                        <bs type="valid" val="904" />                       Bitscore (from blast)
                        <no type="valid" val="0" />                         Number of “other” non-ACGT bases
                        <start type="info" val="99556" />                   Start location in the contig
                        <stop type="info" val="100057" />                   Stop location in the contig
                        <rs type="valid" val="0" />                         Repeat score: sequence identity of secondary hit
                        <al type="valid" val="501" />                       Alignment lengh (from blast)
                        <fwd type="info" val="1" />                         Fwd or rev
                        <cid type="info" val="Contig_34_43.018_pilon" />    Contig identifier
                        <nm type="valid" val="0" />                         Number of mismatches
                        <ngo type="valid" val="0" />                        Number of gaps open
                        <rss type="valid" val="N" />                        Requires start/stop
                        <bc type="valid" val="N" />                         Begin codon (=start codon)
                        <ec type="valid" val="N" />                         End codon
                        <fl type="valid" val="Y" />                         Full-length alignment
                        <is type="valid" val="N" />                         Internal stop
                        <s>ATTGC</s>                                        Sequence
                        <sb>CTGCGTGA</sb>                                   Before the sequence
                        <sa>GCCATCTC</sa>                                   After the sequence
                    </a>
                </as>
            </l>
            ...
        </loci>

        The parser iterate over the different tags, <l>, then <as>, then <a> and finally all the tags in <a>.
        """
        root = fromstring(xml_content)
        for locus in root.iter("l"):
            locus_id = locus.find("id").text
            allele_calls: list[AlleleCall] = []
            for allele in locus.iter("a"):
                data: dict[str, str] = {}
                for node in allele:
                    # type can be "info" or "valid", or absent (for sequences)
                    node_type = node.attrib.get("type")
                    if node_type in ("info", "valid"):
                        data[node.tag] = node.attrib["val"]
                    elif node.tag != "id":
                        data[node.tag] = node.text.strip().upper()
                allele_calls.append(
                    AlleleCall(
                        seq=data["s"],
                        seq_up=data.get("sb", ""),
                        seq_down=data.get("sa", ""),
                        identity=data["si"],
                        evalue=data["ev"],
                        bitscore=data["bs"],
                        align_len=data["al"],
                        has_complete_alignment=data["fl"] == "Y",
                        repeat_score=data["rs"],
                        mismatch_count=data["nm"],
                        gap_count=data["ngo"],
                        contig=data["cid"],
                        pos_start=data["start"],
                        pos_end=data["stop"],
                        is_strand_forward=data["fwd"] == "1",
                        non_atgc_count=data["no"],
                        has_codon_start=data["bc"] == "Y",
                        has_codon_stop=data["ec"] == "Y",
                        is_start_stop_codon_required=data["rss"] == "Y",
                        has_internal_codon_stop=data["is"] == "Y",
                    )
                )

            self.calls[locus_id] = allele_calls

    def get_calls(self) -> dict[str, list[AlleleCall]]:
        return self.calls


class SAMdata:
    def __init__(self, sorted_bam_filename: str, logger: Logger):
        self.sorted_bam_filename = sorted_bam_filename
        self.logger = logger

        self.stats = defaultdict(int)
        # Make sure all the Flag values are present at minimum
        for flag in Flag:
            self.stats[flag.value]

        self.group_to_loci_map: dict[str, set[str]] = {
            LocusGroups.CORE: set(),
            LocusGroups.ACCESSORY: set(),
        }
        self.loci_by_pos: dict[str, dict[int, list[tuple[str, int]]]] = defaultdict(
            lambda: defaultdict(list)
        )

    def bam_to_change(
        self,
        bam_sorted_in: FilePath,
    ):
        self.bam_sorted_in = bam_sorted_in

        self.logger.info("Loading previous BAM calls")

        pysam_index(
            str(self.bam_sorted_in), str(self.bam_sorted_in.with_suffix(".bai"))
        )

        with AlignmentFile(self.bam_sorted_in, "rb") as reader:
            contigs = reader.references
            for record in reader.fetch(until_eof=True):
                locus = record.query_name

                combined_flag = CombinedFlag(record.flag)

                tags_to_values = {
                    tag[0]: list(tag[1:])
                    for tag in record.get_tags(with_value_type=True)
                }
                values = tags_to_values.get("RK", [])
                ranking = values[0] if values else 0

                group = (
                    LocusGroups.CORE
                    if combined_flag.is_flag_set(Flag.CORE)
                    else LocusGroups.ACCESSORY
                )
                self.group_to_loci_map[group].add(locus)

                if ranking != 0:
                    contig = contigs[record.reference_id]
                    positions = record.get_reference_positions()
                    for pos in positions:
                        self.loci_by_pos[contig][pos].append((locus, ranking))

    def change_bam(
        self,
        calls_with_flags: dict[str, dict[int, list[Flag]]],
    ):
        with AlignmentFile(str(self.bam_sorted_in), "rb") as reader:
            with AlignmentFile(
                self.sorted_bam_filename, "wb", template=reader
            ) as writer:
                for record in reader.fetch(until_eof=True):
                    new_tags_counter: dict[str, int] = {
                        flag.tag: 0
                        for flag in (
                            Flag.COVERAGE_PROBLEMS,
                            Flag.BIAS_PROBLEMS,
                            Flag.DOUBLE_CALLS,
                            Flag.MORE_CALLS,
                            Flag.MISMATCH_CALLS,
                        )
                    }
                    tags = record.get_tags(with_value_type=True)
                    tags_to_values = {tag: values for tag, *values in tags}

                    values = tags_to_values.get("RK", [])
                    ranking = values[0] if values else 0
                    if ranking != 0:
                        locus = record.query_name

                        ranking_data = calls_with_flags[locus][ranking]

                        flags_set: set[Flag] = set()
                        for flag in ranking_data:
                            new_tags_counter[flag.tag] += 1
                            flags_set.add(flag)

                        if flags_set:
                            combined_flag = CombinedFlag(record.flag)
                            combined_flag.unset_flags(Flag.CALLED)
                            flags_set.add(Flag.NOT_CALLED)
                            combined_flag.set_flags(*flags_set)
                            record.flag = combined_flag.value

                    tags += [
                        (tag, count, "i") for tag, count in new_tags_counter.items()
                    ]
                    record.set_tags(tags)

                    writer.write(record)
                    self.__save_stats(record.flag)

    def calls_to_bam(
        self,
        calls: dict[str, list[AlleleCall]],
        assembly: FastaInput,
        min_similarity: float,
        core_loci: list[str] = [],
    ):
        self.logger.info(("Calling alleles and writting BAM calls"))

        contigs_seq = get_contigs_seq(str(assembly))
        contigs = sorted(contigs_seq.keys())

        header = {
            "HD": {"VN": "1.0"},  # SAM format version
            "SQ": [
                {"SN": contig, "LN": len(contigs_seq[contig])} for contig in contigs
            ],
        }

        with NamedTemporaryFile() as unsorted_bam_IO:  # Must still exists when calling pysam
            with AlignmentFile(unsorted_bam_IO.name, "wbu", header=header) as writer:
                for locus_id, allele_calls in calls.items():
                    locus_group_flag = (
                        Flag.CORE if locus_id in core_loci else Flag.ACCESSORY
                    )

                    records: list[AlignedSegment] = []
                    if allele_calls:
                        for allele_rank, allele_call in enumerate(
                            allele_calls, start=1
                        ):
                            combined_flag = CombinedFlag()

                            combined_flag.set_flags(locus_group_flag)

                            if allele_call.identity < min_similarity:
                                combined_flag.set_flags(
                                    Flag.IDENTITY_PROBLEMS, Flag.NOT_CALLED
                                )
                            if allele_call.non_atgc_count > 0:
                                combined_flag.set_flags(
                                    Flag.SEQ_PROBLEMS, Flag.NOT_CALLED
                                )
                            if allele_call.has_codon_issues:
                                combined_flag.set_flags(
                                    Flag.CODON_PROBLEMS, Flag.NOT_CALLED
                                )
                            if allele_call.repeat_score > 0:
                                combined_flag.set_flags(
                                    Flag.REPEAT_PROBLEMS, Flag.NOT_CALLED
                                )

                            if not combined_flag.is_flag_set(Flag.NOT_CALLED):
                                combined_flag.set_flags(Flag.CALLED)

                            record = AlignedSegment()
                            record.query_name = locus_id
                            record.reference_id = contigs.index(allele_call.contig)
                            record.reference_start = allele_call.pos_start
                            record.query_sequence = allele_call.seq
                            record.cigartuples = [(0, len(allele_call.seq))]
                            record.set_tags(
                                [
                                    ("ID", str(allele_call.id), "Z"),
                                    ("FI", str(allele_call.flanking_id), "Z"),
                                    ("FS", allele_call.flanking_seq, "Z"),
                                    ("RK", allele_rank, "i"),
                                ]
                            )
                            record.flag = combined_flag.value

                            records.append(record)
                    else:
                        combined_flag = CombinedFlag()
                        combined_flag.set_flags(Flag.NOT_FOUND, locus_group_flag)

                        record = AlignedSegment()
                        record.query_name = locus_id
                        record.flag = combined_flag.value

                        records.append(record)

                    for record in records:
                        writer.write(record)
                        self.__save_stats(record.flag)

            pysam_sort(unsorted_bam_IO.name, "-o", str(self.sorted_bam_filename))

    def __save_stats(self, combined_flag_value: int):
        for flag in Flag:
            self.stats[flag.value] += (combined_flag_value & flag.value) != 0
        self.stats[combined_flag_value] += 1

    def get_formatted_stats(self) -> list[dict[str, str | int | list[dict[str, str]]]]:
        stats = []
        for flag, counts in self.stats.items():
            flag: int
            stats.append(
                {
                    "flag": flag,
                    "counts": counts,
                    "descriptions": [
                        {
                            "tag": Flag.get_flag(1 << i).name,
                            "description": Flag.get_flag(1 << i).description,
                        }
                        for i in range((flag).bit_length())
                        if flag & (1 << i)  # Decompose flag into binary constituants
                    ],
                }
            )
        return stats

    def get_calls(
        self, all_flags: list[Flag] = [], any_flags: list[Flag] = []
    ) -> dict[str, list[FlaggedAlleleCall]]:
        if all_flags and any_flags:
            raise ValueError("Either use all_flags or any_flags")

        with NamedTemporaryFile() as filtered_bam_IO:
            if all_flags or any_flags:
                filtered_bam_filename = filtered_bam_IO.name

                flags = all_flags if all_flags else any_flags
                pysam_view(
                    str(self.sorted_bam_filename),
                    "-o",
                    filtered_bam_filename,
                    "-b",
                    "-u",
                    "--no-header",
                    "--require-flags" if all_flags else "--include-flags",
                    str(sum(flag.value for flag in flags)),
                    catch_stdout=False,
                )
            else:
                filtered_bam_filename = str(self.sorted_bam_filename)

            pysam_index(filtered_bam_filename, filtered_bam_filename + ".bai")

            tag_cov = Flag.COVERAGE_PROBLEMS.tag
            tag_bias = Flag.BIAS_PROBLEMS.tag
            tag_double = Flag.DOUBLE_CALLS.tag
            tag_more = Flag.MORE_CALLS.tag
            tag_mismatch = Flag.MISMATCH_CALLS.tag

            calls: dict[str, list[dict]] = {}
            with AlignmentFile(filtered_bam_filename, "rb") as reader:
                for record in reader.fetch(until_eof=False):
                    if record.query_name not in calls:
                        calls[record.query_name] = []
                    tags = {tag: tag_value for tag, tag_value in record.get_tags()}
                    calls[record.query_name].append(
                        FlaggedAlleleCall(
                            id=int(tags["ID"]),
                            seq=record.query_sequence,
                            flag=record.flag,
                            flank_id=int(tags["FI"]),
                            flank_seq=tags["FS"],
                            coverage_prob=tags.get(tag_cov, None),
                            bias_prob=tags.get(tag_bias, None),
                            double_calls=tags.get(tag_double, None),
                            more_calls=tags.get(tag_more, None),
                            mismatch_calls=tags.get(tag_mismatch, None),
                        )
                    )

        return calls


class Pileup:
    def __init__(
        self,
        cram: FilePath,
        assembly: FastaInput,
        logger: Logger,
        compute_depth_problems_flag: bool = False,
        depth_min: int = 3,
        strand_depth_min: float = 0.33,
        single_nt_call_min: float = 0.7,
        double_nt_call_min: float = 0.8,
        do_double_nt_calling: bool = True,
        use_depth_total: bool = False,
        n_threads: int = 1,
    ) -> None:
        self.assembly = assembly
        self.compute_depth_problems_flag = compute_depth_problems_flag
        self.depth_min = depth_min
        self.strand_depth_min = strand_depth_min
        self.single_nt_call_min = single_nt_call_min
        self.double_nt_call_min = double_nt_call_min
        self.do_double_nt_calling = do_double_nt_calling
        self.use_depth_total = use_depth_total
        self.n_threads = n_threads
        self.logger = logger

        self.logger.info(("Converting CRAM to BAM"))
        self.bam_sorted_in = NamedTemporaryFile()
        pysam_view(
            str(cram),
            "-o",
            self.bam_sorted_in.name,
            "-b",
            "-u",
            "-T",
            str(self.assembly),
            catch_stdout=False,
        )

        pysam_index(self.bam_sorted_in.name, self.bam_sorted_in.name + ".bai")

        self.contigs_seq: dict[str, str] = get_contigs_seq(assembly=self.assembly)
        self.counts: dict[str, dict[str, np.array[np.int32]]] = {
            contig: {} for contig in self.contigs_seq
        }
        self.pos_with_flags: dict[str, dict[int, list[Flag]]] = {}
        self.multiple_calls_flags_counts: dict[str, dict[str, int]] = {}

        msg = "Computing pileups"
        if self.compute_depth_problems_flag:
            msg += " and getting positions with depth problems flags"
        self.logger.info(msg)

        self.__run_parallel_pileup()

        self.bam_sorted_in.close()

    def get_counts(self) -> dict[str, dict[str, np.array]]:
        return self.counts

    def get_pos_with_flags(self) -> dict[str, dict[int, list[Flag]]]:
        return self.pos_with_flags

    def get_multiple_calls_flags_counts(self) -> dict[str, dict[str, int]]:
        return self.multiple_calls_flags_counts

    def get_contigs_len(self) -> dict[str, int]:
        return {contig: len(seq) for contig, seq in self.contigs_seq.items()}

    def get_assembly_len(self) -> int:
        return sum(self.get_contigs_len().values())

    def __run_parallel_pileup(self):
        inputs = [
            {
                "contig": contig,
                "seq": seq,
                "compute_depth_problems_flag": self.compute_depth_problems_flag,
                "depth_min": self.depth_min,
                "strand_depth_min": self.strand_depth_min,
                "single_nt_call_min": self.single_nt_call_min,
                "double_nt_call_min": self.double_nt_call_min,
                "do_double_nt_calling": self.do_double_nt_calling,
                "use_depth_total": self.use_depth_total,
            }
            for contig, seq in sorted(
                self.contigs_seq.items(), key=lambda item: len(item[1]), reverse=True
            )
        ]

        compute_pileup_per_contig_with_files = partial(
            pileup_per_contig,
            bam_sorted_in=self.bam_sorted_in.name,
            assembly=self.assembly,
        )

        results = run_parallel(
            func=compute_pileup_per_contig_with_files,
            inputs=inputs,
            n_threads=self.n_threads,
            logger=self.logger,
        )
        x_mpileup = []
        for (
            contig,
            counts,
            _x_mpileup,
            pos_with_flags,
            multiple_calls_flags_counts,
        ) in results:
            self.counts[contig] = counts
            self.pos_with_flags[contig] = pos_with_flags
            self.multiple_calls_flags_counts[contig] = multiple_calls_flags_counts
            x_mpileup.extend(_x_mpileup)

        if len(x_mpileup) > 0:
            raise ValueError(
                f"Found unrecognized characters in reads pileup: {Counter(x_mpileup)}"
            )


def pileup_per_contig(
    contig: str,
    seq: str,
    bam_sorted_in: str,
    assembly: Path,
    logger: Logger,  # required because passed as keyword during the run_parallel
    compute_depth_problems_flag: bool = False,
    depth_min: int = 3,
    strand_depth_min: float = 0.33,
    single_nt_call_min: float = 0.7,
    double_nt_call_min: float = 0.8,
    do_double_nt_calling: bool = True,
    use_depth_total: bool = False,
):
    def parse_reads(
        pileupreads: list[PileupRead], save_intron_refskip: bool = False
    ) -> dict[str, dict | list]:
        nts = BASIC_NTS.union({"N", "del"})
        if save_intron_refskip:
            nts.union({"refskip"})

        base_counts: dict[str, int] = {nt: 0 for nt in nts}
        counts: dict[str, base_counts | list] = {
            strand.lower(): dict(base_counts) for strand in STRANDS
        }
        counts["X"] = []

        for pileupread in pileupreads:
            alignment = pileupread.alignment

            if pileupread.is_refskip:
                if save_intron_refskip:
                    nt = "refskip"
                else:
                    continue
            elif pileupread.is_del:
                nt = "del"
            else:
                nt = alignment.query_sequence[pileupread.query_position]

            if nt not in nts:
                counts["X"].append(nt)
            else:
                if alignment.is_reverse:
                    strand = "rev"
                else:
                    strand = "fwd"

                counts[strand][nt] += 1

        return counts

    counts: dict[str, np.array[np.int32]] = {
        f"{strand}{nt}": np.zeros(len(seq), dtype=np.int32)
        for nt in BASIC_NTS.union({"Total"})
        for strand in STRANDS
    }

    pos_with_flags: dict[int, list[Flag]] = {}
    multiple_calls_flags_counts: dict[str, int] = {
        flag.name: 0 for flag in MULTIPLE_CALLS_FLAGS
    }
    x_mpileup = []
    with AlignmentFile(bam_sorted_in, "rb", reference_filename=str(assembly)) as reader:
        for pileupcolumn in reader.pileup(
            contig=contig,
            stepper="nofilter",
            flag_filter=1028,
            min_base_quality=13,
            adjust_capq_threshold=0,
            compute_baq=False,
        ):
            pileupcolumn.set_min_base_quality(0)
            pos = pileupcolumn.reference_pos
            _counts = parse_reads(
                pileupreads=pileupcolumn.pileups, save_intron_refskip=False
            )
            for strand in STRANDS:
                for nt in BASIC_NTS:
                    counts[f"{strand}{nt}"][pos] = _counts[f"{strand.lower()}"][f"{nt}"]
            counts["FwdTotal"][pos] = sum(_counts["fwd"].values())
            counts["RevTotal"][pos] = sum(_counts["rev"].values())
            x_mpileup.extend(_counts["X"])

            if compute_depth_problems_flag:
                flags = get_depth_problems_flag(
                    counts=_counts,
                    reference_nt=seq[pos],
                    depth_min=depth_min,
                    strand_depth_min=strand_depth_min,
                    single_nt_call_min=single_nt_call_min,
                    double_nt_call_min=double_nt_call_min,
                    do_double_nt_calling=do_double_nt_calling,
                    use_depth_total=use_depth_total,
                )
                if len(flags) > 0:
                    pos_with_flags[pos] = flags
                    for flag in MULTIPLE_CALLS_FLAGS:
                        if flag in flags:
                            multiple_calls_flags_counts[flag.name] += 1

    return contig, counts, x_mpileup, pos_with_flags, multiple_calls_flags_counts


class FilterWithReads:
    def __init__(
        self,
        cram: Path,
        assembly: Path,
        logger: Logger,
        depth_min: int = 3,
        strand_depth_min: float = 33,
        single_nt_call_min: float = 70,
        double_nt_call_min: float = 80,
        do_double_nt_calling: bool = True,
        use_depth_total: bool = False,
        n_threads: int = 1,
    ) -> None:
        if depth_min < 2:
            raise ValueError(f"depth_min cannot be lower than 2. Provided: {depth_min}")
        if 0 >= strand_depth_min > 50:
            raise ValueError(
                f"strand_depth_min has to be between ]0, 50]. Provided: {strand_depth_min}"
            )
        if 0 >= single_nt_call_min > 100:
            raise ValueError(
                f"single_nt_call_min has to be between ]0, 100]. Provided: {single_nt_call_min}"
            )
        if do_double_nt_calling and single_nt_call_min >= double_nt_call_min > 100:
            raise ValueError(
                f"double_nt_call_min has to be between ]{single_nt_call_min}, 100]. Provided: {double_nt_call_min}"
            )

        strand_depth_min = strand_depth_min / 100
        single_nt_call_min = single_nt_call_min / 100
        double_nt_call_min = double_nt_call_min / 100

        pileup = Pileup(
            cram=cram,
            assembly=assembly,
            n_threads=n_threads,
            logger=logger,
            compute_depth_problems_flag=True,
            depth_min=depth_min,
            strand_depth_min=strand_depth_min,
            single_nt_call_min=single_nt_call_min,
            double_nt_call_min=double_nt_call_min,
            do_double_nt_calling=do_double_nt_calling,
            use_depth_total=use_depth_total,
        )

        logger.info(("Getting positions with coverage problems Flags"))
        self.pos_with_flags = pileup.get_pos_with_flags()
        self.multiple_calls_flags_counts = {
            flag.name: 0 for flag in MULTIPLE_CALLS_FLAGS
        }
        for (
            multiple_calls_flags_counts
        ) in pileup.get_multiple_calls_flags_counts().values():
            for flag_name, flag_count in multiple_calls_flags_counts.items():
                self.multiple_calls_flags_counts[flag_name] += flag_count
        self.assembly_len = pileup.get_assembly_len()

    def get_assembly_len(self) -> int:
        return self.assembly_len

    def get_multiple_calls_flags_counts(self) -> dict[str, int]:
        return self.multiple_calls_flags_counts

    def get_loci_flag(
        self, loci_by_pos: dict[str, dict[int, list[tuple[str, int]]]]
    ) -> dict[str, dict[int, list[Flag]]]:
        calls_with_flags: dict[str, dict[int, list[Flag]]] = defaultdict(
            lambda: defaultdict(list)
        )

        for contig, contig_loci_data in loci_by_pos.items():
            for pos, pos_loci_data in contig_loci_data.items():
                for locus, ranking in pos_loci_data:
                    flags = self.pos_with_flags.get(contig, {}).get(pos)
                    if flags is not None:
                        calls_with_flags[locus][ranking].extend(flags)

        return calls_with_flags


def get_depth_problems_flag(
    counts: dict[str, dict[str, int]],
    reference_nt: str,
    depth_min: int = 3,
    strand_depth_min: float = 0.33,
    single_nt_call_min: float = 0.7,
    double_nt_call_min: float | None = None,
    do_double_nt_calling: bool = True, ## added to fix allele filter error
    use_depth_total: bool = False,
) -> list[Flag]:
    if use_depth_total:
        return calculate_flag_using_total_depth(
            counts,
            reference_nt,
            depth_min,
            strand_depth_min,
            single_nt_call_min,
            double_nt_call_min,
            do_double_nt_calling, ## added to fix allele filter error
        )
    return calculate_flag_per_nt(
        counts,
        reference_nt,
        depth_min,
        strand_depth_min,
        single_nt_call_min,
        double_nt_call_min,
        do_double_nt_calling, ## added to fix allele filter error
    )


def calculate_flag_using_total_depth(
    counts: dict[str, dict[str, int]],
    reference_nt: str,
    depth_min: int,
    strand_depth_min: float,
    single_nt_call_min: float,
    double_nt_call_min: float | None,
    do_double_nt_calling: bool, ## added to fix allele filter error
) -> list[Flag]:
    sums_of_basic_nt_counts = [
        sum(counts[strand.lower()][nt] for nt in BASIC_NTS) for strand in STRANDS
    ]
    sum_of_all_strands = sum(sums_of_basic_nt_counts)
    if sum_of_all_strands < depth_min:
        return [Flag.COVERAGE_PROBLEMS]
    if any(
        s < max(floor(sum_of_all_strands * strand_depth_min), 1)
        for s in sums_of_basic_nt_counts
    ):
        return [Flag.BIAS_PROBLEMS]
    return evaluate_possible_nt(
        reference_nt,
        {nt: sum(counts[strand.lower()][nt] for strand in STRANDS) for nt in BASIC_NTS},
        single_nt_call_min,
        double_nt_call_min,
        do_double_nt_calling, ## added to fix allele filter error
    )


def calculate_flag_per_nt(
    counts: dict[str, dict[str, int]],
    reference_nt: str,
    depth_min: int,
    strand_depth_min: float,
    single_nt_call_min: float,
    double_nt_call_min: float | None,
    do_double_nt_calling: bool, ## added to fix allele filter error
) -> list[Flag]:
    possible_nt_total = {
        nt: [counts[strand.lower()][nt] for strand in STRANDS]
        for nt in BASIC_NTS
        if sum(counts[strand.lower()][nt] for strand in STRANDS) >= depth_min
    }
    if not possible_nt_total:
        return [Flag.COVERAGE_PROBLEMS]

    possible_nt = {
        nt: sum(nt_data)
        for nt, nt_data in possible_nt_total.items()
        if all(s >= max(floor(sum(nt_data) * strand_depth_min), 1) for s in nt_data)
    }
    if not possible_nt:
        return [Flag.BIAS_PROBLEMS]

    return evaluate_possible_nt(
        reference_nt, possible_nt, single_nt_call_min, double_nt_call_min, do_double_nt_calling
    )


def evaluate_possible_nt(
    reference_nt: str,
    possible_nt: dict[str, int],
    single_nt_call_min: float,
    double_nt_call_min: float,
    do_double_nt_calling: bool,
) -> list[Flag]:
    flags: list[Flag] = []

    if len(possible_nt) < 1:
        return flags

    max_nt_count = max(possible_nt.values())
    nt_max = [nt for nt, count in possible_nt.items() if count == max_nt_count]
    if len(nt_max) > 1 or nt_max[0] != reference_nt:
        flags.append(Flag.MISMATCH_CALLS)

    sorted_nt_counts = sorted(possible_nt.values(), reverse=True)
    possible_nt_total = sum(sorted_nt_counts)
    if sorted_nt_counts[0] / possible_nt_total < single_nt_call_min:
        if (
            do_double_nt_calling
            and sum(sorted_nt_counts[:2]) / possible_nt_total >= double_nt_call_min
        ):
            flags.append(Flag.DOUBLE_CALLS)
        else:
            flags.append(Flag.MORE_CALLS)

    return flags
