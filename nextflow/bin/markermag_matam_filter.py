#!/usr/bin/env python3
"""
Wrapper for the MATAM --filter_only step (16S read extraction).

Passes forward and reverse reads to MATAM as separate files, runs
matam_assembly.py --filter_only, then splits MATAM's interleaved aligned
paired FASTQ into predictable R1/R2 outputs.

Called by the MATAM_FILTER Nextflow process.
"""

import os
import sys
import glob
import gzip
import argparse
import subprocess


def sequence_suffix(path):
    name = os.fspath(path).lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    raise ValueError(f"MATAM paired reads must be FASTQ: {path}")


def find_filtered_interleaved(workdir):
    """Find MATAM's paired SortMeRNA FASTX output.

    MATAM passes an explicit ``--aligned`` prefix to SortMeRNA.  With
    ``--paired_in`` and no ``--out2``, SortMeRNA v7 writes one interleaved
    FASTQ at that prefix.  The ``.sortmerna_vs_`` segment is part of MATAM's
    own output naming scheme.
    """
    patterns = (
        "**/*.sortmerna_vs_*.fastq",
        "**/*.sortmerna_vs_*.fastq.gz",
        "**/*.sortmerna_vs_*.fq",
        "**/*.sortmerna_vs_*.fq.gz",
    )
    matches = sorted({
        path
        for pattern in patterns
        for path in glob.glob(os.path.join(workdir, pattern), recursive=True)
        if os.path.isfile(path)
    })
    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one MATAM paired aligned FASTQ under "
            f"{workdir}; found {matches}. MATAM must use SortMeRNA "
            "--paired_in with its standard --aligned output prefix.")
    return matches[0]


def open_fastq(path, mode):
    path = os.fspath(path)
    return gzip.open(path, mode) if path.lower().endswith(".gz") else open(path, mode)


def read_fastq_record(handle, source):
    """Read one FASTQ record, rejecting truncated or malformed input."""
    record = [handle.readline() for _ in range(4)]
    if not record[0]:
        if any(record[1:]):
            raise ValueError(f"Truncated FASTQ record in {source}")
        return None
    if any(not line for line in record):
        raise ValueError(f"Truncated FASTQ record in {source}")
    if not record[0].startswith("@") or not record[2].startswith("+"):
        raise ValueError(f"Malformed FASTQ record in {source}")
    return record


def split_interleaved_fastq(input_path, forward_output, reverse_output):
    """Split an interleaved, pair-preserving FASTQ into R1 and R2 files."""
    pairs = 0
    with open_fastq(input_path, "rt") as source, \
            open_fastq(forward_output, "wt") as forward, \
            open_fastq(reverse_output, "wt") as reverse:
        while True:
            forward_record = read_fastq_record(source, input_path)
            if forward_record is None:
                return pairs
            reverse_record = read_fastq_record(source, input_path)
            if reverse_record is None:
                raise ValueError(
                    f"Interleaved FASTQ has an unpaired final record: {input_path}")
            forward.writelines(forward_record)
            reverse.writelines(reverse_record)
            pairs += 1


def run_filter(args):
    sequence_suffix(args.r1)
    sequence_suffix(args.r2)
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

    filtered_interleaved = find_filtered_interleaved(wd)
    suffix = sequence_suffix(filtered_interleaved)
    forward_output = f"{args.prefix}_16S_R1{suffix}"
    reverse_output = f"{args.prefix}_16S_R2{suffix}"
    split_interleaved_fastq(
        filtered_interleaved, forward_output, reverse_output)
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
