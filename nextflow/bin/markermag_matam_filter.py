#!/usr/bin/env python3
"""
Wrapper for the MATAM --filter_only step (16S read extraction).

Merges R1/R2 with seqtk mergepe, runs matam_assembly.py --filter_only,
then moves the extracted 16S reads file to a predictable output name.

Called by the MATAM_FILTER Nextflow process.
"""

import os
import sys
import glob
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Run MATAM filter-only step for 16S read extraction')
    parser.add_argument('--r1',       required=True,
                        help='Forward reads (fasta or fastq)')
    parser.add_argument('--r2',       required=True,
                        help='Reverse reads (fasta or fastq)')
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
    parser.add_argument('--seqtk',    default='seqtk',
                        help='Path to seqtk (default: seqtk)')
    args = parser.parse_args()

    input_ext = os.path.splitext(args.r1)[1]
    merged    = f'{args.prefix}_merged{input_ext}'
    wd        = f'{args.prefix}_filter_wd'

    # ── Step 1: interleave paired reads for MATAM ─────────────────────────────
    rc = os.system(f'{args.seqtk} mergepe {args.r1} {args.r2} > {merged}')
    if rc != 0:
        print(f'ERROR: seqtk mergepe failed (exit {rc})', file=sys.stderr)
        sys.exit(1)

    # ── Step 2: MATAM --filter_only ───────────────────────────────────────────
    rc = os.system(
        f'{args.matam} -i {merged} -o {wd} '
        f'--cpu {args.cpu} --max_memory {args.mem_mb} -v --filter_only '
        f'-d {args.matam_db}'
    )
    os.remove(merged)

    if rc != 0:
        print(f'ERROR: MATAM --filter_only failed (exit {rc})', file=sys.stderr)
        sys.exit(1)

    # ── Step 3: locate extracted 16S reads inside MATAM workdir ───────────────
    pattern = f'{wd}/workdir/*{input_ext}'
    matches = glob.glob(pattern)
    if not matches:
        print(f'ERROR: No 16S reads found matching {pattern}', file=sys.stderr)
        sys.exit(1)

    out_file = f'{args.prefix}_16S_reads{input_ext}'
    os.rename(matches[0], out_file)
    print(f'16S reads written to: {out_file}')


if __name__ == '__main__':
    main()
