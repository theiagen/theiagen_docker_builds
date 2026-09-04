from ngs_pipeline_lib.base.file import CsvFile, JsonFile, TextFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class AlleleCommonOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self.bam_calls_loci = TextFile(
            name="allele_calls", extension=".bam", compress=False
        )
        self.calls_pcr = JsonFile(name="allele_calls")
        self.stats_calls = JsonFile(name="stats_calls")
        self.dict_standard = JsonFile(name="calls_standard")
        self.csv_core_standard = CsvFile(name="calls_core_standard")
        self.csv_core_pcr = CsvFile(name="calls_core_pcr")
        self.csv_accessory_standard = CsvFile(name="calls_accessory_standard")
        self.csv_accessory_pcr = CsvFile(name="calls_accessory_pcr")


class AlleleCallingOutputs(AlleleCommonOutputs):
    def __init__(self):
        super().__init__()
        self.allele_calls = TextFile(name="allele_calls", extension=".xml")


class AlleleFilteringOutputs(AlleleCommonOutputs):
    def __init__(self):
        super().__init__()
