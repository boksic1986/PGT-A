#!/biosoftware/miniconda/envs/snakemake_env/bin/python
from __future__ import annotations

import argparse


def parse_args():
    from pgta.qc.batch_qc import add_gc_correction_args, add_threshold_args

    parser = argparse.ArgumentParser(
        description="Batch baseline BAM QC for single-bin or multiscale leave-one-out uniformity analysis."
    )
    parser.add_argument("--bams", nargs="+", required=True, help="All baseline BAM paths")
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA path for per-bin GC calculation")
    bin_group = parser.add_mutually_exclusive_group(required=True)
    bin_group.add_argument("--bin-size", type=int, help="Single bin size in bases")
    bin_group.add_argument("--bin-sizes", nargs="+", type=int, help="Target bin sizes in bases")
    parser.add_argument("--mapq", type=int, default=30, help="Minimum MAPQ")
    parser.add_argument("--threads", type=int, default=1, help="Worker count used to count BAMs")
    parser.add_argument("--outdir", required=True, help="Output root directory")
    parser.add_argument("--log", default="", help="Optional log file")
    add_threshold_args(parser)
    add_gc_correction_args(parser)
    return parser.parse_args()


def main():
    args = parse_args()
    from pgta.qc.batch_qc import (
        build_gc_correction_config_from_args,
        build_thresholds_from_args,
        run_batch_bam_qc,
    )
    from pgta.qc.sample_qc import setup_logger

    setup_logger(args.log)
    target_bin_sizes = args.bin_sizes if args.bin_sizes else [args.bin_size]
    run_batch_bam_qc(
        bam_paths=args.bams,
        reference_fasta=args.reference_fasta,
        target_bin_sizes=target_bin_sizes,
        mapq=args.mapq,
        threads=args.threads,
        outdir=args.outdir,
        thresholds=build_thresholds_from_args(args),
        gc_config=build_gc_correction_config_from_args(args),
        nest_binsize_subdir=args.bin_sizes is not None,
    )


if __name__ == "__main__":
    main()
