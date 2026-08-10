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

    def test_single_end_input_is_not_treated_as_interleaved(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            reads = directory / "single.fastq.gz"
            output = directory / "subset.fastq.gz"
            write_fastq(reads, ["read.1", "read.2", "unpaired"])

            total, written = ASSEMBLE.subsample_single_reads(
                str(reads), str(output), 100)

            self.assertEqual((total, written), (3, 3))
            self.assertEqual(
                read_ids(output), ["read.1", "read.2", "unpaired"])


class FilterOutputTests(unittest.TestCase):
    def test_finds_sortmerna_out2_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "workdir" / "sortmerna" / "out"
            output.mkdir(parents=True)
            forward = output / "aligned_fwd.fastq.gz"
            reverse = output / "aligned_rev.fastq.gz"
            write_fastq(forward, ["f1"])
            write_fastq(reverse, ["r1"])

            result = FILTER.find_filtered_pair(directory)

            self.assertEqual(result, (str(forward), str(reverse)))

    def test_rejects_single_combined_filter_output(self):
        with tempfile.TemporaryDirectory() as directory:
            combined = Path(directory) / "aligned.fastq"
            write_fastq(combined, ["f1", "r1"])

            with self.assertRaisesRegex(RuntimeError, "--out2 --paired_in"):
                FILTER.find_filtered_pair(directory)

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
                output = (
                    directory / "sample_filter_wd" / "workdir" /
                    "sortmerna_reads_mapping" / "out")
                output.mkdir(parents=True)
                write_fastq(
                    output / "aligned_paired_fwd.fastq.gz", ["f1"])
                write_fastq(
                    output / "aligned_paired_rev.fastq.gz", ["r1"])

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
                    "MarkerMAG.subsample_reads",
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
