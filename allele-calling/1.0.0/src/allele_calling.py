from pathlib import Path
from shutil import move

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.tools.tools import gunzip_file

from src.inputs import AlleleCallingInputs
from src.outputs import AlleleCallingOutputs
from src.results import ResultsMixin
from src.tools import LocusGroups, parse_reference_loci
from src.transformers import SAMdata, XMLdata


class AlleleCalling(Algorithm[AlleleCallingInputs, AlleleCallingOutputs], ResultsMixin):
    outputs_class = AlleleCallingOutputs

    def execute_stub(self):
        # TODO
        pass

    def execute_implementation(self):
        """
        Calls BLASTAlleleFinder and then parses its output
        """
        self.logger.info(("Prepare input files"))

        working_dir = Path("work")
        working_dir.mkdir(parents=True, exist_ok=True)

        assembly = gunzip_file(self.inputs.assembly, (working_dir / "assembly.fasta"))

        loci = self.get_loci()

        self.run_allele_calling(assembly=assembly, working_dir=working_dir)

        sam_data = self.get_all_calls(
            loci=loci,
            assembly=assembly,
        )

        self.write_calls_output(group_to_loci_map=loci, sam_data=sam_data)

        metrics = self.get_metrics(group_to_loci_map=loci, sam_data=sam_data)
        self.apply_quality_control(metrics)

    def run_allele_calling(self, assembly: Path, working_dir: Path) -> bool:
        (working_dir / "logs").mkdir(exist_ok=True)
        (working_dir / "results").mkdir(exist_ok=True)

        cmd = ["BLASTAlleleFinder_internal"] ## Container - requires full path when not using container
        cmd += ["-_debug_"]
        cmd += ["-_noencrypt_"]
        cmd += ["--nThreads", str(self.inputs.n_threads)]
        cmd += ["--k", "11"]
        cmd += ["--ungapped", "0"]
        cmd += ["--resultsdir", str(working_dir.resolve())]
        cmd += ["--tempdir", str(working_dir.resolve())]
        cmd += ["--localdir", str(working_dir.resolve())]
        cmd += ["--percentIdentity", str(self.inputs.blast_kb.similarity)]
        cmd += ["--query", str(assembly.resolve())]
        cmd += ["--db", str((self.inputs.blast_kb.db / "alleles").resolve())]
        cmd += [
            "--alleleInfo",
            str((self.inputs.blast_kb.db / "alleleinfo.txt").resolve()),
        ]

        accepted_alleles_link_file: Path = (
            self.inputs.blast_kb.db / "acceptedalleles_link"
        )
        if accepted_alleles_link_file.is_file():
            with open(
                accepted_alleles_link_file, encoding="utf-8"
            ) as accepted_alleles_link_reader:
                accepted_alleles_link = accepted_alleles_link_reader.read(1)
        else:
            accepted_alleles_link = "0"
        allele_cache_file = (
            self.inputs.blast_kb.path
            / f"acceptedalleles_{accepted_alleles_link}.fasta.gz"
        )
        if allele_cache_file.is_file():
            self.logger.info("Using accepted alleles.")

        process_stdout, _ = run_external(
            command_line=cmd,
            logger=self.logger,
            working_dir=working_dir,
            use_mamba_env=True, ## looks for executable in container here: "/opt/conda/condabin/" - change back to True for containerization
        )
        self.logger.info(process_stdout)

        move(working_dir / "results" / "loci.xml", self.outputs.allele_calls.path)

    def get_loci(self) -> dict[LocusGroups, set[str]]:
        """
        Return core & accessory locus ids
        """
        with open(self.inputs.blast_kb.loci, encoding="utf-8") as reader:
            lines = []
            # Ignore first line
            next(reader)
            for line in reader:
                lines.append(line)
        core_ids, accessory_ids = parse_reference_loci(lines)
        return {LocusGroups.CORE: core_ids, LocusGroups.ACCESSORY: accessory_ids}

    def get_all_calls(
        self, assembly: Path, loci: dict[LocusGroups, list[str]]
    ) -> SAMdata:
        """
        Extract "normal" and PCR calls
        """

        xml_calls = XMLdata(self.outputs.allele_calls.path).get_calls()

        sam_data = SAMdata(
            sorted_bam_filename=self.outputs.bam_calls_loci.path, logger=self.logger
        )
        sam_data.calls_to_bam(
            calls=xml_calls,
            assembly=assembly,
            min_similarity=self.inputs.blast_kb.similarity,
            core_loci=loci[LocusGroups.CORE],
        )

        return sam_data
