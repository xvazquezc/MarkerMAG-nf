"""Utilities for paired reads that do not depend on renamed identifiers."""

import gzip
import re

from Bio import SeqIO


_MATE_SUFFIX = re.compile(r"/(?:1|2)$")


def canonical_read_id(read_id):
    """Return the pair identifier shared by conventional R1/R2 names."""
    return _MATE_SUFFIX.sub("", read_id.split()[0])


def sam_read_id_and_mate(fields):
    """Return a canonical read ID and mate number from a SAM record.

    Mate identity comes exclusively from SAM FLAG bits 0x40 and 0x80.
    """
    if len(fields) < 2:
        raise ValueError("Malformed SAM record: fewer than two fields")
    flag = int(fields[1])
    is_r1 = bool(flag & 0x40)
    is_r2 = bool(flag & 0x80)
    if is_r1 == is_r2:
        raise ValueError(
            f"SAM record {fields[0]!r} does not identify exactly one mate "
            "with FLAG 0x40/0x80")
    return canonical_read_id(fields[0]), "1" if is_r1 else "2"


def infer_sequence_format(path):
    name = str(path).lower()
    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"
    if name.endswith((".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")):
        return "fasta"
    raise ValueError(f"Unsupported sequence filename: {path}")


def open_sequences(path, mode="rt"):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def extract_reads_by_pair_id(input_path, pair_ids, output_path, sequence_format=None):
    """Stream records whose canonical pair ID occurs in *pair_ids*."""
    wanted = set(pair_ids)
    input_format = sequence_format or infer_sequence_format(input_path)
    output_format = "fastq" if input_format == "fastq" else "fasta-2line"
    written = 0
    with open_sequences(input_path) as input_handle, \
            open_sequences(output_path, "wt") as output_handle:
        for record in SeqIO.parse(input_handle, input_format):
            if canonical_read_id(record.id) in wanted:
                SeqIO.write(record, output_handle, output_format)
                written += 1
    return written


def _read_fastq_record(handle, source):
    """Read one 4-line FASTQ record; None at end of file."""
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


def copy_fastq_pair(forward_input, reverse_input, forward_output, reverse_output):
    """Copy an R1/R2 FASTQ pair, checking that both files stay in step.

    Mates are matched by position and must share a canonical read ID, so
    both CASAVA (identical IDs) and legacy /1, /2 naming are accepted.
    Returns the number of pairs copied.
    """
    pairs = 0
    with open_sequences(forward_input) as forward_source, \
            open_sequences(reverse_input) as reverse_source, \
            open_sequences(forward_output, "wt") as forward, \
            open_sequences(reverse_output, "wt") as reverse:
        while True:
            forward_record = _read_fastq_record(forward_source, forward_input)
            reverse_record = _read_fastq_record(reverse_source, reverse_input)
            if forward_record is None and reverse_record is None:
                return pairs
            if forward_record is None or reverse_record is None:
                raise ValueError(
                    "Paired FASTQ files have different record counts: "
                    f"{forward_input}, {reverse_input}")
            if canonical_read_id(forward_record[0][1:]) != canonical_read_id(reverse_record[0][1:]):
                raise ValueError(
                    "Paired FASTQ records are out of step: "
                    f"{forward_record[0].strip()} vs {reverse_record[0].strip()}")
            forward.writelines(forward_record)
            reverse.writelines(reverse_record)
            pairs += 1
