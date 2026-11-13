#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import pandas as pd

from seqsender_ops.files import (
    combine_passed,
    stage_samples,
    upload_reads,
    write_tarball,
)
from seqsender_ops.metadata import (
    build_sample_metadata,
    write_sample_metadata,
)
from seqsender_ops.models import (
    DATABASES,
    BuildOptions,
    SampleError,
)
from seqsender_ops.validation import (
    run_prep,
    validate_sample_metadata,
    write_metadata_template,
    write_validation_summary,
)
from seqsender_ops.logging_config import setup_logging

logger = logging.getLogger("build_and_prep")


def parse_args():
    parser = argparse.ArgumentParser(description="Builds, prepares, and validates submission files for seqsender")
    parser.add_argument("--metadata_file", type=Path, required=True, help="Terra table TSV, named <table_name>.tsv")
    parser.add_argument("--config_file", type=Path, required=True, help="Path to the seqsender configuration file")
    parser.add_argument("--sample_ids", nargs="+", required=True, help="Sample IDs to be submitted")
    parser.add_argument("--bioproject", type=str, required=True, help="NCBI BioProject ID to submit to")

    # Terra lists these in the same order as the sample IDs
    parser.add_argument("--assembly_fastas", nargs="*", type=Path, default=[], help="Required for GenBank; one assembly FASTA per sample")
    parser.add_argument("--read1_fastqs", nargs="*", type=Path, default=[], help="Required for SRA; one 'read1' FASTQ per sample")
    parser.add_argument("--read2_fastqs", nargs="*", type=Path, default=[], help="Optional for SRA; leave out for single-end reads; one 'read2' FASTQ per sample")

    # with a transfer bucket, NCBI fetches the reads from the bucket ('cloud'); without one, seqsender uploads them over FTP ('local')
    parser.add_argument("--sra_transfer_bucket", type=str, help="Optional; NCBI-readable bucket the passing samples' reads are copied to")

    parser.add_argument("--seqsender_organism_tag", type=str, required=True, help="seqsender organism; options include 'FLU', 'COV', or 'OTHER'")
    parser.add_argument("--biosample", action="store_true")
    parser.add_argument("--sra", action="store_true")
    parser.add_argument("--genbank", action="store_true")
    parser.add_argument("--validation_only", action="store_true", help="Validate only; don't copy the reads to --sra_transfer_bucket")
    parser.add_argument("--output_dir", type=Path, default=Path.cwd(), help="Directory to write the submission files to; populates 'builds' and 'submit' directories")
    return parser.parse_args()


def check_inputs(args) -> None:
    """Exits on problems with the inputs as a whole; problems with a single sample only fail that sample"""
    if not any(getattr(args, db.name) for db in DATABASES):
        raise SystemExit("ERROR: choose at least one of --biosample, --sra or --genbank")
    if args.sra and not args.read1_fastqs:
        raise SystemExit("ERROR: --read1_fastqs is required for SRA submission")
    if args.genbank and not args.assembly_fastas:
        raise SystemExit("ERROR: --assembly_fastas is required for GenBank submission")
    if args.genbank and args.seqsender_organism_tag not in ("FLU", "COV"):
        raise SystemExit(f"ERROR: GenBank submission is only supported for --seqsender_organism_tag FLU and COV, not {args.seqsender_organism_tag}")
    if args.sra_transfer_bucket and not args.sra_transfer_bucket.startswith("gs://"):
        raise SystemExit(f"ERROR: --sra_transfer_bucket must be a gs:// bucket, not {args.sra_transfer_bucket}")

    for name in ("assembly_fastas", "read1_fastqs", "read2_fastqs"):
        files = getattr(args, name)
        if files and len(files) != len(args.sample_ids):
            raise SystemExit(f"ERROR: {len(files)} {name} given for {len(args.sample_ids)} sample_ids; they must be the same length and order")

    # read1/2 FASTQ files must be unique so they don't get overwritten in 'builds/raw_reads'
    if args.sra:
        names = [read.name for read in args.read1_fastqs + args.read2_fastqs]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise SystemExit(f"ERROR: read files must have unique names; repeated: {', '.join(duplicates)}")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output_dir / "seqsender_build_and_prep.log")
    check_inputs(args)
    options = BuildOptions.from_args(args)
    logger.info(f"Preparing {len(args.sample_ids)} samples for {', '.join(db.name for db in options.target_dbs)}")

    builds_dir = args.output_dir / "builds"
    submit_dir = args.output_dir / "submit"

    # 1. stage only the files the selected databases use
    samples = stage_samples(
        builds_dir=builds_dir,
        sample_ids=args.sample_ids,
        assembly_fastas=args.assembly_fastas if options.genbank else [],
        read1_fastqs=args.read1_fastqs if options.sra else [],
        read2_fastqs=args.read2_fastqs if options.sra else [],
    )

    # 2. the Terra table is the base of each sample's metadata.csv; entity:{table_name}_id holds the sample IDs
    metadata_df = pd.read_csv(args.metadata_file, sep="\t", dtype=str, keep_default_na=False)
    sample_column = f"entity:{args.metadata_file.stem}_id"

    # 3. fill in the columns build_and_prep sets itself; the samples are checked in step 5
    for sample in samples:
        sample.metadata = build_sample_metadata(
            sample=sample,
            metadata_df=metadata_df,
            sample_column=sample_column,
            options=options,
        )

    # 4. the template of required columns is written even if every sample fails validation below
    write_metadata_template(
        metadata_df=pd.concat([sample.metadata for sample in samples]),
        sample_column=sample_column,
        options=options,
        template_file=args.output_dir / f"{args.metadata_file.stem}_template.tsv",
    )

    # 5. validate each sample on its own; a sample that fails is left out of the submission
    for sample in samples:
        try:
            # build_and_prep's own checks, then seqsender's
            validate_sample_metadata(
                sample_metadata_df=sample.metadata,
                options=options,
            )
            write_sample_metadata(
                sample=sample,
                options=options,
            )
            run_prep(
                sample=sample,
                options=options,
            )
            sample.passed = True
        except SampleError as error:
            sample.reason = str(error)

        if sample.passed:
            logger.info(f"{sample.sample_id}: PASS")
        else:
            logger.warning(f"{sample.sample_id}: FAIL: {sample.reason}")

    write_validation_summary(samples, args.output_dir)

    # 6. copy valid samples to the `--sra_transfer_bucket` if appropriate. exit if no valid samples
    valid_samples = [sample for sample in samples if sample.passed]
    if not valid_samples:
        logger.warning("No samples passed seqsender prep; no submission tarball written")
        return

    if options.sra and options.sra_transfer_bucket:
        if options.validation_only:
            logger.info(f"--validation_only: the reads were not copied to {options.sra_transfer_bucket}")
        else:
            upload_reads(
                valid_samples=valid_samples,
                read_bucket=options.sra_transfer_bucket,
            )

    # 7. combine the passing samples into submit/ and tar it for the submit task
    combine_passed(
        valid_samples=valid_samples,
        submit_dir=submit_dir,
        config_file=options.config_file,
        include_reads=not options.sra_transfer_bucket,
    )

    write_tarball(
        submit_dir=submit_dir,
        tarball=args.output_dir / "seqsender_post_validation_checkpoint.tar.gz",
    )


if __name__ == "__main__":
    main()
