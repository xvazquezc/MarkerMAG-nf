#!/usr/bin/env python3
"""
Wrapper for the MATAM --filter_only step (16S read extraction).

Passes forward and reverse reads to MATAM as separate files, runs
matam_assembly.py --filter_only, and copies the selected read pairs that MATAM
publishes as separate R1/R2 files to predictable outputs, checking that the
two files stay in step.

Called by the MATAM_FILTER Nextflow process.
"""

import os
import re
import sys
import gzip
import argparse
import subprocess


def sequence_suffix(path):
    name = os.fspath(path).lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    raise ValueError(f"MATAM paired reads must be FASTQ: {path}")


# matam_assembly.py publishes the read pairs selected by its mapping step as
# <out_dir>/marker_reads_R1.fq[.gz] and <out_dir>/marker_reads_R2.fq[.gz].
MATAM_SELECTED_PAIR = ("marker_reads_R1", "marker_reads_R2")


def find_filtered_pair(workdir):
    """Return MATAM's selected (R1, R2) read files from its output directory."""
    pair = []
    for name in MATAM_SELECTED_PAIR:
        candidates = [os.path.join(workdir, name + suffix) for suffix in (".fq.gz", ".fq")]
        matches = [path for path in candidates if os.path.isfile(path)]
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected exactly one MATAM {name}.fq[.gz] in {workdir}; found "
                f"{matches}. MATAM must publish its selected paired reads as "
                "marker_reads_R1/R2.")
        pair.append(matches[0])
    return tuple(pair)


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


_MATE_SUFFIX = re.compile(r"/[12]$")


def pair_id(header_line):
    """Pair identifier of a FASTQ header: first token without a /1 or /2 suffix."""
    return _MATE_SUFFIX.sub("", header_line[1:].split()[0])


def copy_fastq_pair(forward_input, reverse_input, forward_output, reverse_output):
    """Copy an R1/R2 FASTQ pair, checking that both files stay in step."""
    pairs = 0
    with open_fastq(forward_input, "rt") as forward_source, \
            open_fastq(reverse_input, "rt") as reverse_source, \
            open_fastq(forward_output, "wt") as forward, \
            open_fastq(reverse_output, "wt") as reverse:
        while True:
            forward_record = read_fastq_record(forward_source, forward_input)
            reverse_record = read_fastq_record(reverse_source, reverse_input)
            if forward_record is None and reverse_record is None:
                return pairs
            if forward_record is None or reverse_record is None:
                raise ValueError(
                    "Paired FASTQ files have different record counts: "
                    f"{forward_input}, {reverse_input}")
            if pair_id(forward_record[0]) != pair_id(reverse_record[0]):
                raise ValueError(
                    "Paired FASTQ records are out of step: "
                    f"{forward_record[0].strip()} vs {reverse_record[0].strip()}")
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
        # No --no_require_index: MATAM must use the database's prebuilt
        # SortMeRNA index, not have every task index the database (hours for SILVA).
        "-d", args.matam_db,
    ]
    subprocess.run(command, check=True)

    filtered_forward, filtered_reverse = find_filtered_pair(wd)
    suffix = sequence_suffix(filtered_forward)
    forward_output = f"{args.prefix}_16S_R1{suffix}"
    reverse_output = f"{args.prefix}_16S_R2{suffix}"
    copy_fastq_pair(
        filtered_forward, filtered_reverse, forward_output, reverse_output)
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
