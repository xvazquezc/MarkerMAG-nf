import gzip
import io
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
sys.path.insert(0, str(REPOSITORY))
PAIRING = importlib.import_module("MarkerMAG.read_pairing")
MATAM_16S = importlib.import_module("MarkerMAG.matam_16s")


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

def write_matam_selected_pair(directory, forward_ids, reverse_ids, suffix=".fq.gz"):
    """Write the selected read pair as matam_assembly.py publishes it."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    forward = directory / f"marker_reads_R1{suffix}"
    reverse = directory / f"marker_reads_R2{suffix}"
    write_fastq(forward, forward_ids)
    write_fastq(reverse, reverse_ids)
    return forward, reverse


CASAVA_IDS = ["A00123:8:FC:1:1101:1000:1", "A00123:8:FC:1:1101:1000:2"]
CASAVA_R1 = [f"{read_id} 1:N:0:ACGT" for read_id in CASAVA_IDS]
CASAVA_R2 = [f"{read_id} 2:N:0:ACGT" for read_id in CASAVA_IDS]


class FilterOutputTests(unittest.TestCase):
    def test_finds_matam_selected_read_pair(self):
        for find in (FILTER.find_filtered_pair, MATAM_16S.find_filtered_pair):
            with self.subTest(find=find.__module__), \
                    tempfile.TemporaryDirectory() as directory:
                forward, reverse = write_matam_selected_pair(directory, ["f1"], ["f1"])

                self.assertEqual(find(directory), (str(forward), str(reverse)))

    def test_rejects_missing_selected_read_pair(self):
        for find in (FILTER.find_filtered_pair, MATAM_16S.find_filtered_pair):
            with self.subTest(find=find.__module__), \
                    tempfile.TemporaryDirectory() as directory:
                # An interleaved file is not a valid MATAM output any more.
                legacy = Path(directory) / "workdir" / "raw_R1.sortmerna_vs_SILVA_b10_m10.fq.gz"
                legacy.parent.mkdir()
                write_fastq(legacy, ["f1", "f1"])

                with self.assertRaisesRegex(RuntimeError, "marker_reads_R1"):
                    find(directory)

    def test_copies_casava_named_pair_in_step(self):
        for copy in (FILTER.copy_fastq_pair, PAIRING.copy_fastq_pair):
            with self.subTest(copy=copy.__module__), \
                    tempfile.TemporaryDirectory() as directory:
                forward, reverse = write_matam_selected_pair(directory, CASAVA_R1, CASAVA_R2)
                out_forward = Path(directory) / "sample_16S_R1.fq.gz"
                out_reverse = Path(directory) / "sample_16S_R2.fq.gz"

                pairs = copy(forward, reverse, out_forward, out_reverse)

                self.assertEqual(pairs, 2)
                self.assertEqual(read_ids(out_forward), CASAVA_IDS)
                self.assertEqual(read_ids(out_reverse), CASAVA_IDS)

    def test_accepts_legacy_mate_suffixes(self):
        for copy in (FILTER.copy_fastq_pair, PAIRING.copy_fastq_pair):
            with self.subTest(copy=copy.__module__), \
                    tempfile.TemporaryDirectory() as directory:
                forward, reverse = write_matam_selected_pair(directory, ["p1/1"], ["p1/2"])
                directory = Path(directory)

                self.assertEqual(copy(forward, reverse, directory / "o1.fq", directory / "o2.fq"), 1)

    def test_rejects_pairs_out_of_step_or_with_different_counts(self):
        cases = (
            (["p1", "p2"], ["p1", "p3"], "out of step"),
            (["p1", "p2"], ["p1"], "different record counts"),
        )
        for copy in (FILTER.copy_fastq_pair, PAIRING.copy_fastq_pair):
            for forward_ids, reverse_ids, message in cases:
                with self.subTest(copy=copy.__module__, message=message), \
                        tempfile.TemporaryDirectory() as directory:
                    forward, reverse = write_matam_selected_pair(
                        directory, forward_ids, reverse_ids)
                    directory = Path(directory)

                    with self.assertRaisesRegex(ValueError, message):
                        copy(forward, reverse, directory / "o1.fq", directory / "o2.fq")

    def test_filter_invokes_matam_with_two_input_files(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "raw_R1.fastq.gz"
            reverse = directory / "raw_R2.fastq.gz"
            write_fastq(forward, CASAVA_R1[:1])
            write_fastq(reverse, CASAVA_R2[:1])
            args = SimpleNamespace(
                matam="/opt/matam/matam_assembly.py",
                r1=str(forward), r2=str(reverse), prefix="sample",
                cpu=8, mem_mb=30000, matam_db="/db/SILVA",
            )

            def fake_matam(command, check):
                self.assertTrue(check)
                self.assertIn("--forward", command)
                self.assertIn("--reverse", command)
                self.assertIn("--filter_only", command)
                self.assertNotIn("--no_require_index", command)  # prebuilt index only
                self.assertNotIn("-i", command)
                write_matam_selected_pair(
                    directory / "sample_filter_wd", CASAVA_R1[:1], CASAVA_R2[:1])

            with working_directory(directory), mock.patch.object(
                    FILTER.subprocess, "run", side_effect=fake_matam):
                outputs = FILTER.run_filter(args)
            self.assertEqual(outputs, ("sample_16S_R1.fq.gz", "sample_16S_R2.fq.gz"))
            self.assertEqual(read_ids(directory / outputs[0]), CASAVA_IDS[:1])
            self.assertEqual(read_ids(directory / outputs[1]), CASAVA_IDS[:1])

    def run_assembly(self, directory, scaffolds=None):
        """Run the per-depth wrapper with a fake MATAM; return (command, output, stderr)."""
        forward = directory / "16S_R1.fq.gz"
        reverse = directory / "16S_R2.fq.gz"
        write_fastq(forward, CASAVA_R1[:1])
        write_fastq(reverse, CASAVA_R2[:1])
        commands = []

        def fake_matam(command, check):
            self.assertTrue(check)
            commands.append(command)
            if scaffolds is not None:
                workdir = directory / "sample_matam_pct100_wd" / "workdir"
                workdir.mkdir(parents=True)
                (workdir / "scaffolds.NR.min_500bp.fa").write_text(scaffolds)

        argv = ["markermag_matam_assemble_pct.py", "--forward", str(forward),
                "--reverse", str(reverse), "--pct", "100", "--prefix", "sample",
                "--matam_db", "/db/SILVA", "--cpu", "8", "--mem_mb", "30000",
                "--matam", "/opt/matam/matam_assembly.py"]
        stderr = io.StringIO()
        with working_directory(directory), \
                mock.patch.object(sys, "argv", argv), \
                mock.patch.object(ASSEMBLE.subprocess, "run", side_effect=fake_matam), \
                mock.patch("sys.stderr", stderr), mock.patch("sys.stdout", io.StringIO()):
            ASSEMBLE.main()
        self.assertEqual(len(commands), 1)
        output = (directory / "sample_subsample_100_scaffolds.fasta").read_text()
        return commands[0], output, stderr.getvalue()

    def test_assembly_runs_matam_up_to_scaffolding_with_prebuilt_indexes_only(self):
        with tempfile.TemporaryDirectory() as directory:
            command, output, _ = self.run_assembly(
                Path(directory), ">1 count=3\nACGT\n")
        # Only the scaffolds are used: no abundance calculation.
        self.assertEqual(command[command.index("--stop_after") + 1], "scaffolding")
        self.assertNotIn("--no_require_index", command)  # prebuilt indexes only
        self.assertEqual(command[command.index("-d") + 1], "/db/SILVA")
        self.assertEqual(output, ">sample_subsample_100_1\nACGT\n")

    def test_assembly_without_scaffolds_is_an_empty_depth(self):
        # With --stop_after scaffolding MATAM ends normally with no scaffold.
        with tempfile.TemporaryDirectory() as directory:
            _, output, stderr = self.run_assembly(Path(directory), "")
        self.assertEqual(output, "")
        self.assertIn("assembled no scaffolds at pct=100", stderr)


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
