#!/usr/bin/env python3
"""
Wrapper for the MATAM --filter_only step (16S read extraction).

Passes forward and reverse reads to MATAM as separate files, runs
matam_assembly.py --filter_only, then exposes SortMeRNA's separate aligned
paired outputs under predictable names.

Called by the MATAM_FILTER Nextflow process.
"""

import os
import sys
import glob
import shutil
import argparse
import subprocess


def sequence_suffix(path):
    name = path.lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    raise ValueError(f"MATAM paired reads must be FASTQ: {path}")


def find_filtered_pair(workdir):
    forward_patterns = (
        "**/*aligned*fwd*.fastq",
        "**/*aligned*fwd*.fastq.gz",
        "**/*aligned*fwd*.fq",
        "**/*aligned*fwd*.fq.gz",
    )
    reverse_patterns = tuple(
        pattern.replace("fwd", "rev")
        for pattern in forward_patterns
    )

    def matches(patterns):
        found = {
            path
            for pattern in patterns
            for path in glob.glob(os.path.join(workdir, pattern), recursive=True)
            if os.path.isfile(path)
        }
        return sorted(found)

    forward_matches = matches(forward_patterns)
    reverse_matches = matches(reverse_patterns)
    if len(forward_matches) != 1 or len(reverse_matches) != 1:
        raise RuntimeError(
            "Expected exactly one separate aligned forward/reverse FASTQ pair "
            f"under {workdir}; found forward={forward_matches}, "
            f"reverse={reverse_matches}. The MATAM SortMeRNA call must use "
            "--out2 --paired_in.")
    return forward_matches[0], reverse_matches[0]


def run_filter(args):
    wd = f"{args.prefix}_filter_wd"
    command = [
        args.matam,
        "--forward", args.r1,
        "--reverse", args.r2,
        "-o", wd,
        "--cpu", str(args.cpu),
        "--max_memory", str(args.mem_mb),
        "-v", "--filter_only",
        "-d", args.matam_db,
    ]
    subprocess.run(command, check=True)

    filtered_forward, filtered_reverse = find_filtered_pair(wd)
    forward_output = f"{args.prefix}_16S_R1{sequence_suffix(filtered_forward)}"
    reverse_output = f"{args.prefix}_16S_R2{sequence_suffix(filtered_reverse)}"
    shutil.copy2(filtered_forward, forward_output)
    shutil.copy2(filtered_reverse, reverse_output)
    return forward_output, reverse_output


def main():
    parser = argparse.ArgumentParser(
        description='Run MATAM filter-only step for 16S read extraction')
    parser.add_argument('--r1',       required=True,
                        help='Forward reads (fastq or fastq.gz)')
    parser.add_argument('--r2',       required=True,
                        help='Reverse reads (fastq or fastq.gz)')
    parser.add_argument('--prefix',   required=True,
                        help='Sample / output prefix')
    parser.add_argument('--matam_db', required=True,
                        help='Full path prefix to MATAM reference DB')
    parser.add_argument('--cpu',      required=True, type=int,
                        help='CPU threads to allocate to MATAM')
    parser.add_argument('--mem_mb',   required=True, type=int,
                        help='Maximum memory for MATAM (MB)')
    parser.add_argument('--matam',    default='matam_assembly.py',
                        help='Path to matam_assembly.py (default: matam_assembly.py)')
    args = parser.parse_args()

    try:
        forward_output, reverse_output = run_filter(args)
    except FileNotFoundError:
        print(f"ERROR: MATAM executable not found: {args.matam}",
              file=sys.stderr)
        sys.exit(127)
    except (subprocess.CalledProcessError, RuntimeError, ValueError) as error:
        print(f"ERROR: MATAM filtering failed: {error}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            sys.exit(error.returncode or 1)
        sys.exit(1)

    print(
        f"16S read pair written to: {forward_output}, {reverse_output}")


if __name__ == '__main__':
    main()
