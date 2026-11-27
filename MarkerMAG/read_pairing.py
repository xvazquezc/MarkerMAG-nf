"""Utilities for paired reads that do not depend on renamed identifiers."""

import gzip
import os
from itertools import chain, islice

from Bio import SeqIO


_MATE_SEPARATORS = "/._:"
# Separator of the mate number in R1/R2 read IDs: "/" (r/1 r/2), "." (r.1 r.2),
# "_" or ":"; None when both mates share the ID (CASAVA 1.8+, SRA). The read
# IDs are never changed: only the pair identifier used to match mates is.
LEGACY_MATE_SEPARATOR = "/"


def is_same_pair(id1, id2):
    """True if R1 and R2 read IDs (first header tokens) name the same pair.

    Mates share the ID (CASAVA 1.8+, SRA), or the IDs differ only by a trailing
    mate number 1/2 after the same separator: r/1 r/2, r.1 r.2 (MarkerMAG's
    naming, fastq-dump --readids), r_1 r_2, r:1 r:2. Same rule as MATAM's
    fastq_utils.is_same_pair, which checks the input of every MATAM run.
    """
    return id1 == id2 or (len(id1) == len(id2) > 2 and id1[:-1] == id2[:-1]
                          and id1[-2] in _MATE_SEPARATORS and id1[-1] == "1" and id2[-1] == "2")


def mate_separator_of(id1, id2):
    """Mate separator of one R1/R2 ID pair (None for a shared ID).

    Raises ValueError when the IDs do not name the same pair (is_same_pair).
    """
    if id1 == id2:
        return None
    if is_same_pair(id1, id2):
        return id1[-2]
    raise ValueError(f"R1/R2 reads are out of step: {id1} vs {id2}")


def _read_ids(path, count):
    """First-token IDs of the first *count* records of a FASTA/FASTQ(.gz) file."""
    with open_sequences(path) as handle:
        try:
            sequence_format = infer_sequence_format(path)
        except ValueError:  # unusual name: go by the first record
            sequence_format = "fastq" if handle.read(1) == "@" else "fasta"
            handle.seek(0)
        if sequence_format == "fastq":
            headers = (line for number, line in enumerate(handle) if number % 4 == 0)
        else:
            headers = (line for line in handle if line.startswith(">"))
        return [line[1:].split()[0] for line in islice(headers, count)]


def detect_mate_separator(forward_path, reverse_path, pairs=1000):
    """How a sample's R1/R2 read IDs name their mates, from its first pairs.

    Returns the mate separator ("/", ".", "_" or ":") or None for shared IDs.
    Every sampled pair must agree; a mix, or mates out of step, is an error.
    """
    forward_ids = _read_ids(forward_path, pairs)
    reverse_ids = _read_ids(reverse_path, pairs)
    if not forward_ids or len(forward_ids) != len(reverse_ids):
        raise ValueError(
            f"Cannot detect the mate naming of {forward_path} and {reverse_path}: "
            f"{len(forward_ids)} vs {len(reverse_ids)} reads in the first {pairs} pairs")
    separators = {mate_separator_of(id1, id2) for id1, id2 in zip(forward_ids, reverse_ids)}
    if len(separators) > 1:
        names = sorted("shared ID" if s is None else f"'{s}1/{s}2'" for s in separators)
        raise ValueError(
            f"Mixed mate naming in {forward_path} and {reverse_path}: {', '.join(names)}")
    return separators.pop()


def describe_mate_separator(mate_separator):
    if mate_separator is None:
        return "mates share the read ID"
    return f"mates named by a {mate_separator}1/{mate_separator}2 suffix"


def canonical_read_id(read_id, mate_separator=LEGACY_MATE_SEPARATOR):
    """Pair identifier of a read ID: the first token, without the mate suffix
    of the sample's naming (detect_mate_separator); unchanged for shared IDs."""
    token = read_id.split()[0]
    if (mate_separator is not None and len(token) > 2 and token[-2] == mate_separator
            and token[-1] in "12"):
        return token[:-2]
    return token


def sam_read_id_and_mate(fields, mate_separator=LEGACY_MATE_SEPARATOR):
    """Return a canonical read ID and mate number from a SAM record.

    Mate identity comes exclusively from SAM FLAG bits 0x40 and 0x80; the
    pair identifier follows the sample's mate naming (canonical_read_id).
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
    return canonical_read_id(fields[0], mate_separator), "1" if is_r1 else "2"


# A placeholder unmapped record (FLAG 77: paired, unmapped, first mate) with a
# read ID no real read has. See with_end_record.
_END_OF_SAM_RECORD = "\t".join(
    ["\x00end-of-sam", "77", "*", "0", "0", "*", "*", "0", "0", "*", "*"]) + "\n"


def with_end_record(sam_lines):
    """The lines of a read-sorted SAM, then one placeholder unmapped record.

    The SAM-parsing workers write a read pair when the next pair starts; the
    placeholder makes them write the last real pair too (it was lost). Being
    unmapped (CIGAR '*'), the placeholder itself adds nothing.
    """
    return chain(sam_lines, [_END_OF_SAM_RECORD])


def split_sam_by_pair(sam_path, out_folder, lines_per_piece, mate_separator,
                      prefix="splitted_sam_"):
    """Split a SAM sorted by read name into pieces of about *lines_per_piece*
    records, cutting only between read pairs, so that every pair is parsed
    by one worker whatever the number of pieces. `split -l` cut pairs in two.

    Pairs are identified as the parsing workers do (canonical_read_id with the
    sample's mate naming). Header lines are left out (the workers skip them).
    Returns the paths of the pieces, in order.
    """
    pieces = []
    piece = None
    lines_in_piece = 0
    previous_pair = None
    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            pair = canonical_read_id(line.split("\t", 1)[0], mate_separator)
            if piece is None or (lines_in_piece >= lines_per_piece and pair != previous_pair):
                if piece is not None:
                    piece.close()
                path = os.path.join(out_folder, "%s%06d" % (prefix, len(pieces)))
                pieces.append(path)
                piece = open(path, "w")
                lines_in_piece = 0
            piece.write(line)
            lines_in_piece += 1
            previous_pair = pair
    if piece is not None:
        piece.close()
    return pieces


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


def extract_reads_by_pair_id(input_path, pair_ids, output_path, sequence_format=None,
                             mate_separator=LEGACY_MATE_SEPARATOR):
    """Stream records whose canonical pair ID occurs in *pair_ids*."""
    wanted = set(pair_ids)
    input_format = sequence_format or infer_sequence_format(input_path)
    output_format = "fastq" if input_format == "fastq" else "fasta-2line"
    written = 0
    with open_sequences(input_path) as input_handle, \
            open_sequences(output_path, "wt") as output_handle:
        for record in SeqIO.parse(input_handle, input_format):
            if canonical_read_id(record.id, mate_separator) in wanted:
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

    Mates are matched by position and must name the same pair (is_same_pair):
    identical IDs (CASAVA, SRA) or a trailing /1 /2, .1 .2, _1 _2 or :1 :2.
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
            if not is_same_pair(forward_record[0][1:].split()[0], reverse_record[0][1:].split()[0]):
                raise ValueError(
                    "Paired FASTQ records are out of step: "
                    f"{forward_record[0].strip()} vs {reverse_record[0].strip()}")
            forward.writelines(forward_record)
            reverse.writelines(reverse_record)
            pairs += 1
