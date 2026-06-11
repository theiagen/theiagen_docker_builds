from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    DirectoryPath,
    FastaInput,
    FilePath,
    KnowledgeBase,
    OrganismInput,
    QCKnowledgeBase,
)
from pydantic import Field


class AlleleCallingKB(KnowledgeBase):
    similarity: float | None = Field(
        description="Similarity threshold used for BLASTAlleleFinder."
    )
    db: DirectoryPath | None = Field(
        description="Path to the BLASTAlleleFinder database folder."
    )
    loci: FilePath | None = Field(
        description="Path to the file with the reference loci."
    )



class AlleleCallingInputs(BaseInputs):
    assembly: FastaInput = Field(description="Path to the assembly file.")
    blast_kb: AlleleCallingKB = Field(
        description="Path to the BLASTAlleleFinder knowledge base folder."
    )
    qc_kb: QCKnowledgeBase = Field(
        description="Path to the quality control knowledge base folder."
    )
    organism: OrganismInput = Field(
        description="Organism, either genus or genus and species."
    )
    subspecies: str | None = Field(
        description="Subspecies of organism. Currently used for Listeria only."
    )


class AlleleFilteringKB(KnowledgeBase):
    depth_min: int = Field(
        description="Minimum nucleotide/total vertical depth (x times) to consider a nucleotide/position for base calling.",
        default=5,
    )
    strand_depth_min: float = Field(
        description="Minimum percentage of nucleotide/total read strand depth to consider a nucleotide/position for base calling.",
        default=1,
    )
    single_nt_call_min: float = Field(
        description="Minimum vertical depth percentage of the most frequent valid nucleotide for single base calling."
        + " If not reached, at least double call is calling.",
        default=65,
    )
    double_nt_call_min: float = Field(
        description="If do_double_nt_calling is True and no single base can be called,"
        + " minimum vertical depth percentage of the two most frequent (valid) nucleotides for double base calling.",
        default=80,
    )
    do_double_nt_calling: bool = Field(
        description="Do double base calling (if False, only single and multiple base calls will be reported).",
        default=True,
    )
    use_depth_total: bool = Field(
        description="Use total depth and not nucleotide depth when analysing depth problems.",
        default=False,
    )


class AlleleFilteringInputs(BaseInputs):
    assembly: FastaInput = Field(description="Path to the assembly file.")
    alignment: FilePath = Field(description="Path to the CRAM alignment file.")
    calls_bam: FilePath = Field(description="Path to the BAM calls file.")
    filtering_kb: AlleleFilteringKB = Field(
        description="Path to Allele Filtering knowledge base folder."
    )
    qc_kb: QCKnowledgeBase = Field(
        description="Path to the quality control knowledge base folder."
    )
    organism: OrganismInput = Field(
        description="Organism, either genus or genus and species."
    )
    subspecies: str | None = Field(
        description="Subspecies of organism. Currently used for Listeria only."
    )
