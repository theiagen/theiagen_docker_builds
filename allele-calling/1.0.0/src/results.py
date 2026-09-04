from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
import json

from src.sam_flags import Flag
from src.tools import LocusGroups, calls_to_csv, calls_to_dict, calls_to_json
from src.transformers import SAMdata


class ResultsMixin:

    """
    For use on a Algorithm implementation
    """

    def apply_quality_control(self, metrics):
        if self.inputs.organism.genus == "LISTERIA" and self.inputs.subspecies is not None: 
            f = open(self.inputs.subspecies)
            data = json.load(f)
            subspecies = data["ANI"]["metrics"]["Subspecies"]
            self.inputs.organism.genus = self.inputs.organism.genus + "_" + subspecies
            quality_control = QualityControl(
                qc_dict=self.inputs.qc_kb.qc.get_dict(),
                organism=self.inputs.organism,
                report=self.qc_report,
            )
            quality_control.report.add_metrics(metrics=metrics)
            quality_control.apply(section_name="AlleleCalling", observations=metrics)
        else: 
            quality_control = QualityControl(
                qc_dict=self.inputs.qc_kb.qc.get_dict(),
                organism=self.inputs.organism,
                report=self.qc_report,
            )
            quality_control.report.add_metrics(metrics=metrics)
            quality_control.apply(section_name="AlleleCalling", observations=metrics)

    @staticmethod
    def get_metrics(
        group_to_loci_map: dict[LocusGroups, list[str]], sam_data: SAMdata
    ) -> dict[str, int | float]:
        metrics: dict = {}

        n_called_core = len(sam_data.get_calls(all_flags=[Flag.CALLED, Flag.CORE]))
        n_called_accessory = len(
            sam_data.get_calls(all_flags=[Flag.CALLED, Flag.ACCESSORY])
        )

        #Check if no core loci in loci list (ex: Vibrio vulnificus)
        core_loci_list = len(group_to_loci_map[LocusGroups.CORE])
        if core_loci_list == 0:
            core_percent = 0 #if no core loci in list, then set to 0
        else:
            core_percent = round(n_called_core / len(group_to_loci_map[LocusGroups.CORE]) * 100, 2)

        accessory_loci_list = len(group_to_loci_map[LocusGroups.ACCESSORY])
        if accessory_loci_list == 0:
            accessory_percent = 0
        else:
            accessory_percent = round(n_called_accessory / len(group_to_loci_map[LocusGroups.ACCESSORY]) * 100, 2)

        metrics = {
            "coreCount": n_called_core,
            "corePercentage": core_percent,
            "accessoryCount": n_called_accessory,
            "accessoryPercentage": accessory_percent,
            "totalLociCount": n_called_core + n_called_accessory
        }

        return metrics

    def write_calls_output(
        self, sam_data: SAMdata, group_to_loci_map: dict[str, set[str]]
    ):
        self.logger.info(("Writing output files"))

        self.outputs.stats_calls.content = {
            "stats_calls": sam_data.get_formatted_stats()
        }

        calls_called = sam_data.get_calls(all_flags=[Flag.CALLED])
        calls_all = sam_data.get_calls()

        core_loci = group_to_loci_map[LocusGroups.CORE]
        acc_loci = group_to_loci_map[LocusGroups.ACCESSORY]
        all_loci = core_loci.union(acc_loci)
        self.outputs.dict_standard.content = calls_to_dict(
            sample_id=self.inputs.sample_id, calls=calls_called, loci=all_loci
        )

        self.outputs.calls_pcr.content = calls_to_json(
            sample_id=self.inputs.sample_id, calls=calls_all, loci=group_to_loci_map
        )

        self.outputs.csv_core_standard.content = calls_to_csv(
            sample_id=self.inputs.sample_id,
            calls=calls_called,
            loci=core_loci,
        )
        self.outputs.csv_core_pcr.content = calls_to_csv(
            sample_id=self.inputs.sample_id,
            calls=calls_all,
            loci=core_loci,
        )

        self.outputs.csv_accessory_standard.content = calls_to_csv(
            sample_id=self.inputs.sample_id,
            calls=calls_called,
            loci=acc_loci,
        )
        self.outputs.csv_accessory_pcr.content = calls_to_csv(
            sample_id=self.inputs.sample_id,
            calls=calls_all,
            loci=acc_loci,
        )
