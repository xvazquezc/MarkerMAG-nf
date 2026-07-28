#!/usr/bin/env python3
"""
Wrapper for one MATAM assembly run at a specific subsample percentage.

Subsamples the extracted 16S reads (using helpers from MarkerMAG.matam_16s),
then runs matam_assembly.py for this percentage only and prefixes the output
scaffolds with the sample/subsample tag.

This script is called once per subsample percentage by the MATAM_ASSEMBLE
Nextflow process, enabling all percentages to run as independent HPC jobs.
"""

import os
import sys
import argparse


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
    parser.add_argument('--usearch',   default='usearch',
                        help='Path to usearch (default: usearch)')
    parser.add_argument('--seqtk',     default='seqtk',
                        help='Path to seqtk (default: seqtk)')
    args = parser.parse_args()

    # Import subsample helpers from the installed MarkerMAG package
    from MarkerMAG.matam_16s import subsample_sortmerna_output, prefix_seq

    # Derive a clean tag (int when pct is a whole number, float otherwise)
    pct_tag    = int(args.pct) if args.pct == int(args.pct) else args.pct
    input_ext  = os.path.splitext(args.reads_16s)[1]
    subsampled = f'{args.prefix}_16S_subset_{pct_tag}{input_ext}'
    matam_wd   = f'{args.prefix}_matam_pct{pct_tag}_wd'
    assemblies = f'{matam_wd}/workdir/scaffolds.NR.min_500bp.fa'
    out_fasta  = f'{args.prefix}_subsample_{pct_tag}_scaffolds.fasta'

    # ── Step 1: subsample 16S reads at this percentage ────────────────────────
    subsample_sortmerna_output(
        args.reads_16s, args.pct, subsampled, args.usearch, args.seqtk)

    # ── Step 2: run MATAM assembly ────────────────────────────────────────────
    rc = os.system(
        f'{args.matam} -d {args.matam_db} -i {subsampled} '
        f'--cpu {args.cpu} --max_memory {args.mem_mb} -v -o {matam_wd}'
    )
    if rc != 0:
        print(
            f'WARNING: MATAM assembly returned non-zero exit for pct={pct_tag}',
            file=sys.stderr)

    # ── Step 3: prefix scaffolds and write named output ───────────────────────
    if os.path.isfile(assemblies):
        prefix_seq(assemblies, f'{args.prefix}_subsample_{pct_tag}', out_fasta)
        print(f'Assemblies written to: {out_fasta}')
    else:
        # No scaffolds produced at this depth — non-fatal; Nextflow output is optional
        print(
            f'WARNING: No scaffolds produced by MATAM at pct={pct_tag}. '
            f'This subsampling depth may be too low.',
            file=sys.stderr)
        sys.exit(0)


if __name__ == '__main__':
    main()
