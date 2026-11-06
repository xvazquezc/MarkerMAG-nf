#!/usr/bin/env python3
"""
Wrapper for one MATAM assembly run at a specific subsample percentage.

Deterministically subsamples two separate paired-read files, then runs
matam_assembly.py for this percentage and prefixes the output scaffolds with
the sample/subsample tag.

This script is called once per subsample percentage by the MATAM_ASSEMBLE
Nextflow process, enabling all percentages to run as independent HPC jobs.
"""

import os
import sys
import gzip
import random
import argparse
import subprocess
from itertools import zip_longest

from Bio import SeqIO


def sequence_format(path):
    name = path.lower()
    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"
    raise ValueError(f"Paired MATAM reads must be FASTQ: {path}")


def sequence_suffix(path):
    name = path.lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    raise ValueError(f"Unsupported sequence filename: {path}")


def open_sequences(path, mode="rt"):
    return gzip.open(path, mode) if path.lower().endswith(".gz") else open(path, mode)


def subsample_paired_reads(forward_path, reverse_path,
                           forward_output, reverse_output, pct, seed=42):
    if not 0 < pct <= 100:
        raise ValueError(f"Subsample percentage must be > 0 and <= 100, got {pct}")

    forward_fmt = sequence_format(forward_path)
    reverse_fmt = sequence_format(reverse_path)
    if forward_fmt != reverse_fmt:
        raise ValueError("Forward and reverse reads must use the same format")

    output_fmt = "fastq" if forward_fmt == "fastq" else "fasta-2line"
    rng = random.Random(seed)
    total_pairs = 0
    written_pairs = 0

    with open_sequences(forward_path) as forward_handle, \
            open_sequences(reverse_path) as reverse_handle, \
            open_sequences(forward_output, "wt") as forward_output_handle, \
            open_sequences(reverse_output, "wt") as reverse_output_handle:
        forward_records = SeqIO.parse(forward_handle, forward_fmt)
        reverse_records = SeqIO.parse(reverse_handle, reverse_fmt)
        for forward_record, reverse_record in zip_longest(
                forward_records, reverse_records):
            if forward_record is None or reverse_record is None:
                raise ValueError(
                    "Forward and reverse files contain different numbers of reads")
            total_pairs += 1
            if pct == 100 or rng.random() < pct / 100:
                SeqIO.write(forward_record, forward_output_handle, output_fmt)
                SeqIO.write(reverse_record, reverse_output_handle, output_fmt)
                written_pairs += 1

    return total_pairs, written_pairs


def main():
    parser = argparse.ArgumentParser(
        description='Run MATAM assembly for one subsample percentage')
    parser.add_argument(
        '--forward', required=True,
        help='Forward extracted 16S reads (fastq or fastq.gz)')
    parser.add_argument(
        '--reverse', required=True,
        help='Reverse extracted 16S reads (fastq or fastq.gz)')
    parser.add_argument('--pct',       required=True, type=float,
                        help='Subsample percentage (e.g. 25)')
    parser.add_argument('--seed',      default=42, type=int,
                        help='Deterministic random seed for paired subsampling '
                             '(default: %(default)s)')
    parser.add_argument('--prefix',    required=True,
                        help='Sample / output prefix')
    parser.add_argument('--matam_db',  required=True,
                        help='Full path prefix to MATAM reference DB')
    parser.add_argument('--cpu',       required=True, type=int,
                        help='CPU threads to allocate to MATAM')
    parser.add_argument('--mem_mb',    required=True, type=int,
                        help='Maximum memory for MATAM (MB)')
    parser.add_argument('--matam',     default='matam_assembly.py',
                        help='Path to matam_assembly.py (default: matam_assembly.py)')
    args = parser.parse_args()

    pct_tag = int(args.pct) if args.pct == int(args.pct) else args.pct
    matam_wd = f"{args.prefix}_matam_pct{pct_tag}_wd"
    assemblies = f"{matam_wd}/workdir/scaffolds.NR.min_500bp.fa"
    out_fasta = f"{args.prefix}_subsample_{pct_tag}_scaffolds.fasta"

    forward_subset = (
        f"{args.prefix}_16S_subset_{pct_tag}_R1"
        f"{sequence_suffix(args.forward)}")
    reverse_subset = (
        f"{args.prefix}_16S_subset_{pct_tag}_R2"
        f"{sequence_suffix(args.reverse)}")
    total_records, written_records = subsample_paired_reads(
        args.forward, args.reverse, forward_subset, reverse_subset,
        args.pct, args.seed)
    input_args = [
        "--forward", forward_subset,
        "--reverse", reverse_subset,
    ]
    unit = "pairs"
    destination = f"{forward_subset} and {reverse_subset}"

    print(
        f"Subsampled {written_records} of {total_records} {unit} "
        f"at {args.pct:g}% with seed {args.seed} into {destination}")

    if written_records == 0:
        open(out_fasta, "w").close()
        print(
            f"WARNING: No reads selected at pct={pct_tag}; "
            "emitting an empty assembly for this depth.",
            file=sys.stderr)
        return

    command = [
        args.matam, "-d", args.matam_db,
        *input_args,
        "--cpu", str(args.cpu),
        "--max_memory", str(args.mem_mb), "-v", "-o", matam_wd,
    ]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        print(f"ERROR: MATAM executable not found: {args.matam}", file=sys.stderr)
        sys.exit(127)
    except subprocess.CalledProcessError as error:
        print(
            f"ERROR: MATAM assembly failed at pct={pct_tag} "
            f"with exit status {error.returncode}",
            file=sys.stderr)
        sys.exit(error.returncode or 1)

    if os.path.isfile(assemblies):
        with open(out_fasta, "w") as output_handle:
            for record in SeqIO.parse(assemblies, "fasta"):
                record.id = f"{args.prefix}_subsample_{pct_tag}_{record.id}"
                record.description = ""
                SeqIO.write(record, output_handle, "fasta-2line")
        print(f"Assemblies written to: {out_fasta}")
    else:
        open(out_fasta, "w").close()
        print(
            f"WARNING: MATAM completed but produced no scaffolds at pct={pct_tag}; "
            "emitting an empty assembly for this depth.",
            file=sys.stderr)


if __name__ == '__main__':
    main()
