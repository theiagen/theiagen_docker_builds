import factory
from factory import post_generation

from src.models import AlleleCall, LocusCall


class AlleleCallFactory(factory.Factory):
    class Meta:
        model = AlleleCall

    allele_id = factory.Faker("md5")
    qualities = {}
    info = {}
    sequences = {}
    ss_codons_problems = factory.Faker("boolean")


class EmptyLocusCallFactory(factory.Factory):
    class Meta:
        model = LocusCall

    locus_id = factory.Faker("md5")
    seq_presence_absence = factory.Faker("boolean")
    ignore_start_stop_codons = factory.Faker("boolean")


class LocusCallSingleAlleleFactory(EmptyLocusCallFactory):
    @post_generation
    def post(obj: LocusCall, create, extracted, **kwargs):
        if "codons_problems" in kwargs:
            ss_codons_problems = kwargs["codons_problems"]
        else:
            ss_codons_problems = False
        if "allele_id" in kwargs:
            obj.add_allele_call(
                AlleleCallFactory.build(
                    allele_id=kwargs["allele_id"], ss_codons_problems=ss_codons_problems
                )
            )
        else:
            obj.add_allele_call(
                AlleleCallFactory.build(ss_codons_problems=ss_codons_problems)
            )


class LocusCallMultiAlleleFactory(EmptyLocusCallFactory):
    @post_generation
    def post(obj: LocusCall, create, extracted, **kwargs):
        if "codons_problems" in kwargs:
            ss_codons_problems = kwargs["codons_problems"]
        else:
            ss_codons_problems = False
        for i in range(5):
            if "allele_id" in kwargs:
                obj.add_allele_call(
                    AlleleCallFactory.build(
                        allele_id=kwargs["allele_id"],
                        ss_codons_problems=ss_codons_problems,
                    )
                )
            else:
                obj.add_allele_call(
                    AlleleCallFactory.build(ss_codons_problems=ss_codons_problems)
                )
