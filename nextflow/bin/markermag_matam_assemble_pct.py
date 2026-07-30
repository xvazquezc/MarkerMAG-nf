#!/usr/bin/env python3
"""
Wrapper for one MATAM assembly run at a specific subsample percentage.

Deterministically subsamples the extracted 16S reads, keeping .1/.2 mates
together, then runs matam_assembly.py for this percentage only and prefixes
the output scaffolds with the sample/subsample tag.

This script is called once per subsample percentage by the MATAM_ASSEMBLE
Nextflow process, enabling all percentages to run as independent HPC jobs.
"""

import os
import sys
import gzip
import random
import argparse
import subprocess

from Bio import SeqIO


def sequence_format(path):
    name = path.lower()
    return "fastq" if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")) else "fasta"


def open_sequences(path):
    return gzip.open(path, "rt") if path.lower().endswith(".gz") else open(path)


def read_group(record_id):
    base, separator, mate = record_id.rpartition(".")
    return base if separator and mate in {"1", "2"} else record_id


def subsample_reads(input_path, output_path, pct):
    if not 0 < pct <= 100:
        raise ValueError(f"Subsample percentage must be > 0 and <= 100, got {pct}")

    fmt = sequence_format(input_path)
    rng = random.Random(1)
    selected_groups = {}
    total_records = 0

    with open_sequences(input_path) as input_handle:
        for record in SeqIO.parse(input_handle, fmt):
            total_records += 1
            group = read_group(record.id)
            if group not in selected_groups:
                selected_groups[group] = pct == 100 or rng.random() < pct / 100

    written_records = 0
    output_fmt = "fastq" if fmt == "fastq" else "fasta-2line"
    with open_sequences(input_path) as input_handle, open(output_path, "w") as output_handle:
        for record in SeqIO.parse(input_handle, fmt):
            if selected_groups[read_group(record.id)]:
                SeqIO.write(record, output_handle, output_fmt)
                written_records += 1
    return total_records, written_records


def main():
    parser = argparse.ArgumentParser(
        description='Run MATAM assembly for one subsample percentage')
    parser.add_argument('--reads_16s', required=True,
                        help='Extracted 16S reads (fasta or fastq)')
    parser.add_argument('--pct',       required=True, type=float,
                        help='Subsample percentage (e.g. 25)')
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
    input_fmt = sequence_format(args.reads_16s)
    input_ext = ".fastq" if input_fmt == "fastq" else ".fasta"
    subsampled = f"{args.prefix}_16S_subset_{pct_tag}{input_ext}"
    matam_wd = f"{args.prefix}_matam_pct{pct_tag}_wd"
    assemblies = f"{matam_wd}/workdir/scaffolds.NR.min_500bp.fa"
    out_fasta = f"{args.prefix}_subsample_{pct_tag}_scaffolds.fasta"

    total_records, written_records = subsample_reads(
        args.reads_16s, subsampled, args.pct)
    print(
        f"Subsampled {written_records} of {total_records} records "
        f"at {args.pct:g}% into {subsampled}")

    if written_records == 0:
        open(out_fasta, "w").close()
        print(
            f"WARNING: No reads selected at pct={pct_tag}; "
            "emitting an empty assembly for this depth.",
            file=sys.stderr)
        return

    command = [
        args.matam, "-d", args.matam_db, "-i", subsampled,
        "--cpu", str(args.cpu), "--max_memory", str(args.mem_mb),
        "-v", "-o", matam_wd,
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
