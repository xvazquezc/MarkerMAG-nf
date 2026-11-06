import gzip
import importlib
import importlib.util
import os
import sys
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from Bio import SeqIO


REPOSITORY = Path(__file__).resolve().parents[2]


def load_script(name):
    path = REPOSITORY / "nextflow" / "bin" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ASSEMBLE = load_script("markermag_matam_assemble_pct")
FILTER = load_script("markermag_matam_filter")
PAIRING = importlib.import_module("MarkerMAG.read_pairing")


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def write_fastq(path, ids):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for record_id in ids:
            handle.write(
                f"@{record_id}\n"
                "ACGTACGT\n"
                "+\n"
                "IIIIIIII\n"
            )


def write_interleaved_fastq(path, forward_ids, reverse_ids):
    if len(forward_ids) != len(reverse_ids):
        raise ValueError("interleaved test input requires complete pairs")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for forward_id, reverse_id in zip(forward_ids, reverse_ids):
            for record_id in (forward_id, reverse_id):
                handle.write(
                    f"@{record_id}\n"
                    "ACGTACGT\n"
                    "+\n"
                    "IIIIIIII\n"
                )


def read_ids(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [record.id for record in SeqIO.parse(handle, "fastq")]


class PairedSubsamplingTests(unittest.TestCase):
    def test_gzip_pair_is_kept_in_separate_files_and_by_position(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "reads_R1.fastq.gz"
            reverse = directory / "reads_R2.fastq.gz"
            forward_output = directory / "subset_R1.fastq.gz"
            reverse_output = directory / "subset_R2.fastq.gz"
            forward_ids = ["forward-A", "forward-B", "forward-C"]
            reverse_ids = ["reverse-X", "reverse-Y", "reverse-Z"]
            write_fastq(forward, forward_ids)
            write_fastq(reverse, reverse_ids)

            total, written = ASSEMBLE.subsample_paired_reads(
                str(forward), str(reverse),
                str(forward_output), str(reverse_output), 100)

            self.assertEqual((total, written), (3, 3))
            self.assertEqual(read_ids(forward_output), forward_ids)
            self.assertEqual(read_ids(reverse_output), reverse_ids)

    def test_seed_reproduces_a_subsample_and_can_select_a_different_one(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "reads_R1.fastq"
            reverse = directory / "reads_R2.fastq"
            forward_ids = [f"f{index}" for index in range(20)]
            reverse_ids = [f"r{index}" for index in range(20)]
            write_fastq(forward, forward_ids)
            write_fastq(reverse, reverse_ids)

            def subsample(seed, tag):
                forward_output = directory / f"{tag}_R1.fastq"
                reverse_output = directory / f"{tag}_R2.fastq"
                ASSEMBLE.subsample_paired_reads(
                    str(forward), str(reverse), str(forward_output),
                    str(reverse_output), 50, seed=seed)
                return read_ids(forward_output), read_ids(reverse_output)

            first = subsample(17, "first")
            repeated = subsample(17, "repeated")
            different = subsample(18, "different")

            self.assertEqual(first, repeated)
            self.assertNotEqual(first, different)
            self.assertEqual(
                [identifier.replace("f", "r", 1) for identifier in first[0]],
                first[1])

    def test_pair_with_different_record_counts_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "reads_R1.fastq"
            reverse = directory / "reads_R2.fastq"
            write_fastq(forward, ["f1", "f2"])
            write_fastq(reverse, ["r1"])

            with self.assertRaisesRegex(
                    ValueError, "different numbers of reads"):
                ASSEMBLE.subsample_paired_reads(
                    str(forward), str(reverse),
                    str(directory / "out_R1.fastq"),
                    str(directory / "out_R2.fastq"), 100)

class FilterOutputTests(unittest.TestCase):
    def test_finds_matam_interleaved_sortmerna_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "workdir"
            output.mkdir(parents=True)
            aligned = output / "raw_R1.sortmerna_vs_SILVA_b10_m10.fastq"
            write_interleaved_fastq(aligned, ["f1"], ["r1"])

            result = FILTER.find_filtered_interleaved(directory)

            self.assertEqual(result, str(aligned))

    def test_rejects_missing_matam_aligned_output(self):
        with tempfile.TemporaryDirectory() as directory:
            unrelated = Path(directory) / "raw_R1.interleaved.fastq"
            write_fastq(unrelated, ["f1", "r1"])

            with self.assertRaisesRegex(RuntimeError, "paired aligned FASTQ"):
                FILTER.find_filtered_interleaved(directory)

    def test_splits_interleaved_output_into_positionally_paired_files(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            aligned = directory / "aligned.fastq.gz"
            forward = directory / "filtered_R1.fastq.gz"
            reverse = directory / "filtered_R2.fastq.gz"
            write_interleaved_fastq(aligned, ["f1", "f2"], ["r1", "r2"])

            pairs = FILTER.split_interleaved_fastq(aligned, forward, reverse)

            self.assertEqual(pairs, 2)
            self.assertEqual(read_ids(forward), ["f1", "f2"])
            self.assertEqual(read_ids(reverse), ["r1", "r2"])

    def test_rejects_unpaired_interleaved_output(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            aligned = directory / "aligned.fastq"
            write_fastq(aligned, ["f1"])

            with self.assertRaisesRegex(ValueError, "unpaired final record"):
                FILTER.split_interleaved_fastq(
                    aligned, directory / "filtered_R1.fastq",
                    directory / "filtered_R2.fastq")

    def test_filter_invokes_matam_with_two_input_files(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "raw_R1.fastq.gz"
            reverse = directory / "raw_R2.fastq.gz"
            write_fastq(forward, ["f1"])
            write_fastq(reverse, ["r1"])
            args = SimpleNamespace(
                matam="/opt/matam/matam_assembly.py",
                r1=str(forward), r2=str(reverse), prefix="sample",
                cpu=8, mem_mb=30000, matam_db="/db/SILVA",
            )

            def fake_matam(command, check):
                self.assertTrue(check)
                self.assertIn("--forward", command)
                self.assertIn("--reverse", command)
                self.assertNotIn("-i", command)
                output = directory / "sample_filter_wd" / "workdir"
                output.mkdir(parents=True)
                write_interleaved_fastq(
                    output / "raw_R1.sortmerna_vs_SILVA_b10_m10.fastq",
                    ["f1"], ["r1"])

            with working_directory(directory), mock.patch.object(
                    FILTER.subprocess, "run", side_effect=fake_matam):
                outputs = FILTER.run_filter(args)
            self.assertEqual(read_ids(directory / outputs[0]), ["f1"])
            self.assertEqual(read_ids(directory / outputs[1]), ["r1"])


class PythonCompatibilityTests(unittest.TestCase):
    def test_all_package_modules_import(self):
        sys.path.insert(0, str(REPOSITORY))
        try:
            for module_name in (
                    "MarkerMAG.MarkerMAG_config",
                    "MarkerMAG.barrnap_16s",
                    "MarkerMAG.get_cp_num",
                    "MarkerMAG.link_16s",
                    "MarkerMAG.matam_16s",
                    "MarkerMAG.polish_16s",
                    "MarkerMAG.uclust_16s"):
                with self.subTest(module=module_name):
                    importlib.import_module(module_name)
        finally:
            sys.path.remove(str(REPOSITORY))

    def test_sam_flags_identify_mates_without_name_suffixes(self):
        r1_id, r1_mate = PAIRING.sam_read_id_and_mate(
            ["A00123:42:H7L2YDSX7:1:1101:1000:1000", "65"])
        r2_id, r2_mate = PAIRING.sam_read_id_and_mate(
            ["A00123:42:H7L2YDSX7:1:1101:1000:1000", "129"])
        self.assertEqual(r1_id, r2_id)
        self.assertEqual((r1_mate, r2_mate), ("1", "2"))

    def test_sam_record_without_one_mate_flag_is_rejected(self):
        for flag in ("0", "192"):
            with self.subTest(flag=flag), self.assertRaisesRegex(
                    ValueError, "FLAG 0x40/0x80"):
                PAIRING.sam_read_id_and_mate(["read", flag])

    def test_streaming_extraction_preserves_original_gzip_read_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "input_R1.fastq.gz"
            output = directory / "selected_R1.fastq.gz"
            write_fastq(source, ["instrument/1", "another/1"])

            written = PAIRING.extract_reads_by_pair_id(
                source, {"instrument"}, output)

            self.assertEqual(written, 1)
            self.assertEqual(read_ids(output), ["instrument/1"])
            self.assertEqual(PAIRING.canonical_read_id("sample.1"), "sample.1")


if __name__ == "__main__":
    unittest.main()
