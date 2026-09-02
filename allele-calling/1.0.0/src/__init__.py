__version__ = "4.0.1"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.allele_calling import AlleleCalling
from src.allele_filtering import AlleleFiltering
from src.inputs import AlleleCallingInputs, AlleleFilteringInputs


@cli.command(name="AlleleCalling")
def allele_calling(inputs: AlleleCallingInputs):
    setup_logging(inputs.logging_dir)
    algorithm = AlleleCalling(inputs)
    algorithm.execute()


@cli.command(name="AlleleFiltering")
def allele_filtering(inputs: AlleleFilteringInputs):
    setup_logging(inputs.logging_dir)
    algorithm = AlleleFiltering(inputs)
    algorithm.execute()
