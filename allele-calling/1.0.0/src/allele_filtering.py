from pathlib import Path

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.tools import gunzip_file

from src.inputs import AlleleFilteringInputs
from src.outputs import AlleleFilteringOutputs
from src.results import ResultsMixin
from src.transformers import FilterWithReads, SAMdata


class AlleleFiltering(
    Algorithm[AlleleFilteringInputs, AlleleFilteringOutputs], ResultsMixin
):
    outputs_class = AlleleFilteringOutputs

    def execute_stub(self):
        # TODO
        pass

    def execute_implementation(self):
        self.logger.info(("Prepare input files"))

        working_dir = Path("work")
        working_dir.mkdir(parents=True, exist_ok=True)

        assembly = gunzip_file(self.inputs.assembly, (working_dir / "assembly.fasta"))

        sam_data, assembly_len, multiple_calls_flags_counts = self.run_allele_filtering(
            assembly=assembly
        )

        self.write_calls_output(
            sam_data=sam_data, group_to_loci_map=sam_data.group_to_loci_map
        )

        metrics = self.get_metrics(
            group_to_loci_map=sam_data.group_to_loci_map, sam_data=sam_data
        )
        self.add_multiple_calls_metrics(
            metrics, multiple_calls_flags_counts, assembly_len
        )
        self.apply_quality_control(metrics)

        assembly.unlink()

    def run_allele_filtering(self, assembly: Path) -> SAMdata:
        sam_data = SAMdata(
            sorted_bam_filename=self.outputs.bam_calls_loci.path, logger=self.logger
        )
        sam_data.bam_to_change(bam_sorted_in=self.inputs.calls_bam)

        filter_with_reads = FilterWithReads(
            cram=self.inputs.alignment,
            assembly=assembly,
            depth_min=self.inputs.filtering_kb.depth_min,
            strand_depth_min=self.inputs.filtering_kb.strand_depth_min,
            single_nt_call_min=self.inputs.filtering_kb.single_nt_call_min,
            double_nt_call_min=self.inputs.filtering_kb.double_nt_call_min,
            do_double_nt_calling=self.inputs.filtering_kb.do_double_nt_calling,
            use_depth_total=self.inputs.filtering_kb.use_depth_total,
            n_threads=self.inputs.n_threads,
            logger=self.logger,
        )
        sam_data.change_bam(
            calls_with_flags=filter_with_reads.get_loci_flag(
                loci_by_pos=sam_data.loci_by_pos
            ),
        )

        return (
            sam_data,
            filter_with_reads.get_assembly_len(),
            filter_with_reads.get_multiple_calls_flags_counts(),
        )

    def add_multiple_calls_metrics(
        self,
        metrics: dict,
        multiple_calls_flags_counts: dict[str:int],
        assembly_len: int,
    ):
        multiple_calls_count = 0
        for flag_name, flag_count in multiple_calls_flags_counts.items():
            metrics[f"nt_{flag_name}_count"] = flag_count
            metrics[f"nt_{flag_name}_percentage"] = round(
                flag_count / assembly_len * 100, 4
            )
            multiple_calls_count += flag_count
        metrics["nt_MULTIPLE_CALLS_count"] = multiple_calls_count
        metrics["nt_MULTIPLE_CALLS_percentage"] = round(
            multiple_calls_count / assembly_len * 100, 4
        )
