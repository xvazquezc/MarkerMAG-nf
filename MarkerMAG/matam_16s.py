#!/usr/bin/env python3

import argparse
import gzip
import os
import random
import shutil
import subprocess
from datetime import datetime
from itertools import zip_longest
from shutil import which

from Bio import SeqIO

from MarkerMAG.read_pairing import copy_fastq_pair


matam_16s_usage = """
=============================== matam_16s examples ===============================

MarkerMAG matam_16s -p Test -r1 R1.fastq.gz -r2 R2.fastq.gz -t 12 -d DB_PREFIX
MarkerMAG matam_16s -p Test -r16s1 extracted_R1.fastq.gz \
    -r16s2 extracted_R2.fastq.gz -t 12 -d DB_PREFIX

Paired reads remain in separate files throughout; interleaved input is unsupported.

=================================================================================
"""


def report_and_log(message, log_file, quiet):
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    with open(log_file, "a", encoding="utf-8") as log_handle:
        log_handle.write(f"{timestamp} {message}\n")
    if not quiet:
        print(f"{timestamp} {message}")


def str_to_num_list(values):
    percentages = []
    for value in (float(item) for item in values.split(",")):
        if not 0 < value <= 100:
            raise ValueError(
                f"Subsample percentages must be > 0 and <= 100, got {value}")
        percentages.append(int(value) if value.is_integer() else value)
    return sorted(percentages)


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
    if path.lower().endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def subsample_pair(forward_path, reverse_path, forward_output, reverse_output,
                   percentage):
    forward_type = sequence_format(forward_path)
    reverse_type = sequence_format(reverse_path)
    if forward_type != reverse_type:
        raise ValueError("Forward and reverse reads must use the same format")

    output_type = "fastq"
    rng = random.Random(1)
    total = 0
    written = 0
    with open_sequences(forward_path) as forward_handle, \
            open_sequences(reverse_path) as reverse_handle, \
            open_sequences(forward_output, "wt") as forward_output_handle, \
            open_sequences(reverse_output, "wt") as reverse_output_handle:
        forward_records = SeqIO.parse(forward_handle, forward_type)
        reverse_records = SeqIO.parse(reverse_handle, reverse_type)
        for forward_record, reverse_record in zip_longest(
                forward_records, reverse_records):
            if forward_record is None or reverse_record is None:
                raise ValueError(
                    "Forward and reverse files contain different numbers of reads")
            total += 1
            if percentage == 100 or rng.random() < percentage / 100:
                SeqIO.write(forward_record, forward_output_handle, output_type)
                SeqIO.write(reverse_record, reverse_output_handle, output_type)
                written += 1
    return total, written


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
                f"MATAM filtering must publish exactly one {name}.fq[.gz] in "
                f"{workdir}; found {matches}")
        pair.append(matches[0])
    return tuple(pair)


def prefix_fasta(input_path, prefix, output_path):
    with open(output_path, "w", encoding="utf-8") as output_handle:
        for record in SeqIO.parse(input_path, "fasta"):
            record.id = f"{prefix}_{record.id}"
            record.description = ""
            SeqIO.write(record, output_handle, "fasta-2line")


def parse_uclust_output(uclust_output, cluster_members):
    clusters = {}
    with open(uclust_output, encoding="utf-8") as input_handle:
        for line in input_handle:
            fields = line.rstrip().split("\t")
            cluster_id = int(fields[1])
            clusters.setdefault(cluster_id, set()).add(fields[8].split()[0])
    with open(cluster_members, "w", encoding="utf-8") as output_handle:
        for cluster_id in sorted(clusters):
            members = ",".join(sorted(clusters[cluster_id]))
            output_handle.write(f"Cluster_{cluster_id}\t{members}\n")


def check_inputs(args):
    raw_pair = args["r1"] is not None or args["r2"] is not None
    extracted_pair = args["r16s1"] is not None or args["r16s2"] is not None

    if raw_pair and not (args["r1"] and args["r2"]):
        raise ValueError("-r1 and -r2 must be provided together")
    if extracted_pair and not (args["r16s1"] and args["r16s2"]):
        raise ValueError("-r16s1 and -r16s2 must be provided together")
    if sum((raw_pair, extracted_pair)) != 1:
        raise ValueError(
            "Provide exactly one paired input mode: -r1/-r2 or -r16s1/-r16s2")

    paths = [
        path for path in (
            args["r1"], args["r2"], args["r16s1"], args["r16s2"])
        if path is not None
    ]
    missing = [path for path in paths if not os.path.isfile(path)]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    for path in paths:
        sequence_format(path)
    return raw_pair, extracted_pair


def matam_16s(args):
    raw_pair, extracted_pair = check_inputs(args)
    if which(args["matam"]) is None and not os.path.isfile(args["matam"]):
        raise FileNotFoundError(f"Executable not found: {args['matam']}")
    if which("vsearch") is None:
        raise FileNotFoundError("Executable not found on PATH: vsearch")

    output_prefix = args["p"]
    workdir = f"{output_prefix}_Matam16S_wd"
    if os.path.isdir(workdir):
        if not args["force"]:
            raise FileExistsError(f"Output folder detected: {workdir}")
        shutil.rmtree(workdir)
    os.mkdir(workdir)

    log_file = os.path.join(workdir, f"{output_prefix}_matam_16s.log")
    if raw_pair:
        filter_dir = os.path.join(workdir, f"{output_prefix}_get_16S_reads_wd")
        command = [
            args["matam"], "--forward", args["r1"], "--reverse", args["r2"],
            "-o", filter_dir, "--cpu", str(args["t"]),
            "--max_memory", str(args["mem"]), "-v", "--filter_only",
            "-d", args["d"],
        ]
        report_and_log("Extracting paired 16S reads with MATAM",
                       log_file, args["quiet"])
        subprocess.run(command, check=True)
        filtered_forward, filtered_reverse = find_filtered_pair(filter_dir)
        reads_forward = os.path.join(
            workdir, f"{output_prefix}_16S_R1{sequence_suffix(filtered_forward)}")
        reads_reverse = os.path.join(
            workdir, f"{output_prefix}_16S_R2{sequence_suffix(filtered_forward)}")
        copy_fastq_pair(
            filtered_forward, filtered_reverse, reads_forward, reads_reverse)
    else:
        reads_forward = args["r16s1"]
        reads_reverse = args["r16s2"]

    assemblies = []
    for percentage in str_to_num_list(args["pct"]):
        percentage_tag = f"{percentage:g}" if isinstance(
            percentage, float) else str(percentage)
        output_dir = os.path.join(
            workdir, f"{output_prefix}_16S_subset_{percentage_tag}_Matam_wd")

        forward_subset = os.path.join(
            workdir,
            f"{output_prefix}_16S_subset_{percentage_tag}_R1"
            f"{sequence_suffix(reads_forward)}")
        reverse_subset = os.path.join(
            workdir,
            f"{output_prefix}_16S_subset_{percentage_tag}_R2"
            f"{sequence_suffix(reads_reverse)}")
        total, written = subsample_pair(
            reads_forward, reads_reverse, forward_subset, reverse_subset,
            percentage)
        input_arguments = [
            "--forward", forward_subset, "--reverse", reverse_subset]

        report_and_log(
            f"Selected {written} of {total} "
            f"pairs at {percentage_tag}%",
            log_file, args["quiet"])
        if written == 0:
            continue

        command = [
            args["matam"], "-d", args["d"], *input_arguments,
            "--cpu", str(args["t"]), "--max_memory", str(args["mem"]),
            "-v", "-o", output_dir,
        ]
        subprocess.run(command, check=True)

        # MATAM's stable output name, not the internal workdir file (whose
        # name depends on MATAM's options).
        assembly = os.path.join(output_dir, "final_assembly.fa")
        if os.path.isfile(assembly):
            prefixed = os.path.join(
                output_dir, "workdir", "final_assembly.prefixed.fa")
            prefix_fasta(
                assembly, f"{output_prefix}_subsample_{percentage_tag}",
                prefixed)
            assemblies.append(prefixed)

    combined = os.path.join(
        workdir, f"{output_prefix}_assembled_16S_unclustered.fasta")
    with open(combined, "w", encoding="utf-8") as output_handle:
        for assembly in assemblies:
            with open(assembly, encoding="utf-8") as input_handle:
                shutil.copyfileobj(input_handle, output_handle)

    clustered = os.path.join(
        workdir, f"{output_prefix}_assembled_16S_uclust_{args['i']}.fasta")
    uclust_table = os.path.join(
        workdir, f"{output_prefix}_assembled_16S_uclust_{args['i']}.uc")
    members = os.path.join(
        workdir, f"{output_prefix}_assembled_16S_uclust_{args['i']}.txt")
    subprocess.run([
        "vsearch", "--cluster_fast", combined, "--id", str(args["i"]),
        "--centroids", clustered, "--uc", uclust_table, "--quiet",
    ], check=True)
    parse_uclust_output(uclust_table, members)
    report_and_log(f"Assemblies written to {clustered}",
                   log_file, args["quiet"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=matam_16s_usage)
    parser.add_argument("-p", required=True, help="output prefix")
    parser.add_argument("-r1", default=None, help="forward raw reads")
    parser.add_argument("-r2", default=None, help="reverse raw reads")
    parser.add_argument("-r16s1", default=None,
                        help="forward extracted 16S reads")
    parser.add_argument("-r16s2", default=None,
                        help="reverse extracted 16S reads")
    parser.add_argument("-pct", default="1,5,10,25,50,75,100")
    parser.add_argument("-d", required=True, help="MATAM reference DB prefix")
    parser.add_argument("-i", type=float, default=0.999)
    parser.add_argument("-t", type=int, default=1)
    parser.add_argument("-mem", type=int, default=10240)
    parser.add_argument("-force", action="store_true")
    parser.add_argument("-quiet", action="store_true")
    parser.add_argument("-matam", default="matam_assembly.py")
    matam_16s(vars(parser.parse_args()))
