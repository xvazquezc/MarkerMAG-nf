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
                    "MarkerMAG.rename_reads",
                    "MarkerMAG.subsample_reads",
                    "MarkerMAG.uclust_16s"):
                with self.subTest(module=module_name):
                    importlib.import_module(module_name)
        finally:
            sys.path.remove(str(REPOSITORY))

    def test_gzip_rename_works_with_multiprocessing_context(self):
        rename_reads = importlib.import_module("MarkerMAG.rename_reads")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "input_R1.fastq.gz"
            reverse = directory / "input_R2.fastq.gz"
            write_fastq(forward, ["old-f1", "old-f2"])
            write_fastq(reverse, ["old-r1", "old-r2"])
            with working_directory(directory):
                rename_reads.rename_reads({
                    "r1": str(forward), "r2": str(reverse),
                    "p": "renamed", "fq": True, "t": 2,
                })
            self.assertEqual(
                read_ids(directory / "renamed_R1.fastq"),
                ["renamed_1.1", "renamed_2.1"])
            self.assertEqual(
                read_ids(directory / "renamed_R2.fastq"),
                ["renamed_1.2", "renamed_2.2"])


if __name__ == "__main__":
    unittest.main()
