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


def write_matam_alignments(directory, text="@HD\tVN:1.0\n"):
    """Write the selected reads' alignments as matam_assembly.py publishes them:
    marker_alignments.sam.gz, a link into its workdir."""
    directory = Path(directory)
    target = directory / "workdir" / "reads.sortmerna_vs_SILVA_b10_m10.sam.gz"
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(target, "wt") as handle:
        handle.write(text)
    link = directory / "marker_alignments.sam.gz"
    link.symlink_to("workdir/" + target.name)
    return link


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

    def test_finds_matam_selected_alignments(self):
        for find in (FILTER.find_filtered_alignments, MATAM_16S.find_filtered_alignments):
            with self.subTest(find=find.__module__), \
                    tempfile.TemporaryDirectory() as directory:
                link = write_matam_alignments(directory)

                self.assertEqual(find(directory), str(link))

    def test_rejects_missing_selected_alignments(self):
        # A MATAM without --read_alignments does not publish them.
        for find in (FILTER.find_filtered_alignments, MATAM_16S.find_filtered_alignments):
            with self.subTest(find=find.__module__), \
                    tempfile.TemporaryDirectory() as directory:
                write_matam_selected_pair(directory, ["f1"], ["f1"])

                with self.assertRaisesRegex(RuntimeError, "marker_alignments"):
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

    def test_accepts_every_mate_suffix_matam_accepts(self):
        # .1/.2 is MarkerMAG's own read naming: it used to stop MATAM_FILTER
        # with "Paired FASTQ records are out of step".
        for copy in (FILTER.copy_fastq_pair, PAIRING.copy_fastq_pair):
            for r1, r2 in (("p1/1", "p1/2"), ("Heavy_1691.1", "Heavy_1691.2"),
                           ("p1_1", "p1_2"), ("p1:1", "p1:2"), ("SRR1.7", "SRR1.7")):
                with self.subTest(copy=copy.__module__, r1=r1), \
                        tempfile.TemporaryDirectory() as directory:
                    forward, reverse = write_matam_selected_pair(directory, [r1], [r2])
                    directory = Path(directory)

                    self.assertEqual(
                        copy(forward, reverse, directory / "o1.fq", directory / "o2.fq"), 1)

    def test_pair_rule_matches_matam(self):
        for module in (FILTER, PAIRING):
            for id1, id2, same in (("r", "r", True), ("r/1", "r/2", True), ("r.1", "r.2", True),
                                   ("r.1", "r/2", False), ("r.2", "r.1", False), ("r1.1", "r2.2", False),
                                   ("a.1", "a.1", True), ("x1", "x2", False)):
                with self.subTest(module=module.__name__, ids=(id1, id2)):
                    self.assertIs(module.is_same_pair(id1, id2), same)

    def test_rejects_pairs_out_of_step_or_with_different_counts(self):
        cases = (
            (["p1", "p2"], ["p1", "p3"], "out of step"),
            (["p1.1"], ["p2.2"], "out of step"),
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

    def filter_args(self, forward, reverse, alignments_only=False):
        return SimpleNamespace(
            matam="/opt/matam/matam_assembly.py",
            r1=str(forward), r2=str(reverse), prefix="sample",
            cpu=8, mem_mb=30000, matam_db="/db/SILVA", alignments_only=alignments_only,
        )

    def test_mapping_extracted_reads_writes_only_their_alignments(self):
        """MATAM_MAP_16S: the inputs may carry the names MATAM_FILTER publishes."""
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "sample_16S_R1.fq.gz"
            reverse = directory / "sample_16S_R2.fq.gz"
            write_fastq(forward, CASAVA_R1)
            write_fastq(reverse, CASAVA_R2)

            def fake_matam(command, check):
                write_matam_selected_pair(
                    directory / "sample_filter_wd", CASAVA_R1[:1], CASAVA_R2[:1])
                write_matam_alignments(directory / "sample_filter_wd")

            with working_directory(directory), mock.patch.object(
                    FILTER.subprocess, "run", side_effect=fake_matam):
                outputs = FILTER.run_filter(self.filter_args(forward, reverse, True))
            self.assertEqual(outputs, (None, None, "sample_16S.sam.gz"))
            self.assertEqual(read_ids(forward), CASAVA_IDS)
            self.assertEqual(read_ids(reverse), CASAVA_IDS)

    def test_filter_refuses_to_write_over_its_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "sample_16S_R1.fq.gz"
            reverse = directory / "raw_R2.fq.gz"
            write_fastq(forward, CASAVA_R1)
            write_fastq(reverse, CASAVA_R2)
            with working_directory(directory), mock.patch.object(
                    FILTER.subprocess, "run") as run, \
                    self.assertRaisesRegex(ValueError, "would overwrite the input"):
                FILTER.run_filter(self.filter_args(forward, reverse))
            run.assert_not_called()
            self.assertEqual(read_ids(forward), CASAVA_IDS)

    def test_filter_invokes_matam_with_two_input_files(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            forward = directory / "raw_R1.fastq.gz"
            reverse = directory / "raw_R2.fastq.gz"
            write_fastq(forward, CASAVA_R1[:1])
            write_fastq(reverse, CASAVA_R2[:1])
            args = self.filter_args(forward, reverse)

            def fake_matam(command, check):
                self.assertTrue(check)
                self.assertIn("--forward", command)
                self.assertIn("--reverse", command)
                self.assertIn("--filter_only", command)
                self.assertNotIn("--no_require_index", command)  # prebuilt index only
                self.assertNotIn("-i", command)
                write_matam_selected_pair(
                    directory / "sample_filter_wd", CASAVA_R1[:1], CASAVA_R2[:1])
                write_matam_alignments(directory / "sample_filter_wd", "@HD\tVN:1.0\nrecord\n")

            with working_directory(directory), mock.patch.object(
                    FILTER.subprocess, "run", side_effect=fake_matam):
                outputs = FILTER.run_filter(args)
            self.assertEqual(outputs, ("sample_16S_R1.fq.gz", "sample_16S_R2.fq.gz",
                                       "sample_16S.sam.gz"))
            self.assertEqual(read_ids(directory / outputs[0]), CASAVA_IDS[:1])
            self.assertEqual(read_ids(directory / outputs[1]), CASAVA_IDS[:1])
            # The alignments file itself, not the link MATAM publishes.
            self.assertFalse((directory / outputs[2]).is_symlink())
            with gzip.open(directory / outputs[2], "rt") as handle:
                self.assertEqual(handle.read(), "@HD\tVN:1.0\nrecord\n")

    def run_assembly(self, directory, scaffolds=None, extra_argv=()):
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
                # As MATAM: the scaffolds in workdir/, under a name that
                # depends on its options, and final_assembly.fa linking to them.
                out_dir = directory / "sample_matam_pct100_wd"
                (out_dir / "workdir").mkdir(parents=True)
                (out_dir / "workdir" / "scaffolds.NR.min_800bp.fa").write_text(scaffolds)
                (out_dir / "final_assembly.fa").symlink_to("workdir/scaffolds.NR.min_800bp.fa")

        argv = ["markermag_matam_assemble_pct.py", "--forward", str(forward),
                "--reverse", str(reverse), "--pct", "100", "--prefix", "sample",
                "--matam_db", "/db/SILVA", "--cpu", "8", "--mem_mb", "30000",
                "--matam", "/opt/matam/matam_assembly.py", *extra_argv]
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
        # MATAM's graph compaction is seeded (it seeds from the clock otherwise).
        self.assertEqual(command[command.index("--seed") + 1], "42")

    def test_assembly_passes_the_matam_seed(self):
        with tempfile.TemporaryDirectory() as directory:
            command, _, _ = self.run_assembly(
                Path(directory), ">1\nACGT\n", ("--seed", "5", "--matam_seed", "7"))
        # --seed is the subsampling seed; MATAM gets --matam_seed.
        self.assertEqual(command[command.index("--seed") + 1], "7")
        self.assertEqual(command.count("--seed"), 1)

    def test_assembly_passes_the_edge_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            command, _, _ = self.run_assembly(Path(directory), ">1\nACGT\n")
        self.assertEqual(command[command.index("-E") + 1], "1")  # MATAM's default
        with tempfile.TemporaryDirectory() as directory:
            command, _, _ = self.run_assembly(
                Path(directory), ">1\nACGT\n", ("--min_overlap_edge", "0"))
        self.assertEqual(command[command.index("-E") + 1], "0")
        self.assertEqual(command.count("-E"), 1)

    def test_assembly_rejects_a_negative_edge_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(SystemExit):
                self.run_assembly(Path(directory), ">1\nACGT\n", ("--min_overlap_edge", "-1"))
        with self.assertRaisesRegex(Exception, "must be >= 0"):
            ASSEMBLE.non_negative_int("-1")
        self.assertEqual(ASSEMBLE.non_negative_int("0"), 0)

    def test_assembly_takes_the_read_mapping_from_the_given_alignments(self):
        with tempfile.TemporaryDirectory() as directory:
            command, _, _ = self.run_assembly(
                Path(directory), ">1\nACGT\n", ("--sam", "sample_16S.sam.gz"))
            self.assertEqual(command[command.index("--read_alignments") + 1],
                             "sample_16S.sam.gz")
        with tempfile.TemporaryDirectory() as directory:
            command, _, _ = self.run_assembly(Path(directory), ">1\nACGT\n")
            self.assertNotIn("--read_alignments", command)

    def test_assembly_without_matam_output_is_an_empty_depth(self):
        with tempfile.TemporaryDirectory() as directory:
            _, output, stderr = self.run_assembly(Path(directory))
        self.assertEqual(output, "")
        self.assertIn("produced no scaffolds at pct=100", stderr)

    def test_assembly_without_scaffolds_is_an_empty_depth(self):
        # With --stop_after scaffolding MATAM ends normally with no scaffold.
        with tempfile.TemporaryDirectory() as directory:
            _, output, stderr = self.run_assembly(Path(directory), "")
        self.assertEqual(output, "")
        self.assertIn("assembled no scaffolds at pct=100", stderr)


def make_matam_db(directory, indexes=("record", "record")):
    """A MATAM db prefix with fasta files and the given clustered/complete indexes:
    "record" (MATAM index record), "stats" (legacy index statistics) or None."""
    prefix = Path(directory) / "db" / "SILVA"
    prefix.parent.mkdir()
    for kind, index in zip(("clustered", "complete"), indexes):
        Path(f"{prefix}.{kind}.fasta").write_text(">r\nACGT\n")
        idx_dir = Path(f"{prefix}.{kind}.smr") / "idx"
        idx_dir.mkdir(parents=True)
        if index == "record":
            (idx_dir / "matam_index.json").write_text("{}")
        elif index == "stats":
            (idx_dir / "123.idx_stats.yaml").write_text("")
    return str(prefix)


class StandaloneMatam16STests(unittest.TestCase):
    def test_unindexed_database_is_rejected_before_running(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = make_matam_db(directory, (None, None))
            with self.assertRaisesRegex(RuntimeError, "index_ref_db.py -i " + prefix):
                MATAM_16S.check_matam_db_indexed(prefix)

    def test_each_index_is_checked_and_legacy_indexes_are_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = make_matam_db(directory, ("stats", None))
            with self.assertRaises(RuntimeError) as error:
                MATAM_16S.check_matam_db_indexed(prefix)
            self.assertIn("complete.smr", str(error.exception))
            self.assertNotIn("clustered.smr", str(error.exception))
        with tempfile.TemporaryDirectory() as directory:
            MATAM_16S.check_matam_db_indexed(make_matam_db(directory, ("stats", "record")))

    def run_standalone(self, directory, **extra_args):
        """Run matam_16s on extracted reads with fake tools; return MATAM commands and output."""
        forward, reverse = directory / "16S_R1.fq", directory / "16S_R2.fq"
        write_fastq(forward, CASAVA_R1)
        write_fastq(reverse, CASAVA_R2)
        args = {"p": "s", "r1": None, "r2": None, "r16s1": str(forward),
                "r16s2": str(reverse), "pct": "50,100", "d": make_matam_db(directory),
                "i": 0.999, "t": 2, "mem": 1000, "force": False, "quiet": True,
                "matam": "matam_assembly.py", "subsample_seed": 1, "matam_seed": 7,
                **extra_args}
        matam_commands = []

        def fake_run(command, check):
            if command[0] == "vsearch":
                uc = command[command.index("--uc") + 1]
                Path(uc).write_text("S\t0\t4\t*\t*\t*\t*\t*\tx\t*\n")
                Path(command[command.index("--centroids") + 1]).write_text(">x\nACGT\n")
                return
            matam_commands.append(command)
            out_dir = Path(command[command.index("-o") + 1])
            if "--filter_only" in command:
                write_matam_alignments(out_dir)
                return
            (out_dir / "workdir").mkdir(parents=True)
            (out_dir / "workdir" / "scaffolds.NR.min_800bp.fa").write_text(">1\nACGT\n")
            (out_dir / "final_assembly.fa").symlink_to("workdir/scaffolds.NR.min_800bp.fa")

        with working_directory(directory), \
                mock.patch.object(MATAM_16S, "which", return_value="/bin/true"), \
                mock.patch.object(MATAM_16S.subprocess, "run", side_effect=fake_run):
            MATAM_16S.matam_16s(args)
            combined = Path("s_Matam16S_wd/s_assembled_16S_unclustered.fasta").read_text()
        return matam_commands, combined

    def test_extracted_reads_are_mapped_once_for_all_subsamples(self):
        with tempfile.TemporaryDirectory() as directory:
            commands, _ = self.run_standalone(Path(directory))
        mapping, *assemblies = commands
        self.assertIn("--filter_only", mapping)
        self.assertEqual(len(assemblies), 2)
        alignments = "s_Matam16S_wd/s_map_16S_reads_wd/marker_alignments.sam.gz"
        for command in assemblies:
            self.assertEqual(command[command.index("--read_alignments") + 1], alignments)

    def test_map_each_subsample_maps_every_subsample_again(self):
        with tempfile.TemporaryDirectory() as directory:
            commands, _ = self.run_standalone(Path(directory), map_each_subsample=True)
        self.assertEqual(len(commands), 2)
        for command in commands:
            self.assertNotIn("--filter_only", command)
            self.assertNotIn("--read_alignments", command)

    def test_assembly_runs_stop_after_scaffolding_seeded_and_read_final_assembly(self):
        with tempfile.TemporaryDirectory() as directory:
            commands, combined = self.run_standalone(Path(directory))
        assemblies = [command for command in commands if "--filter_only" not in command]
        self.assertTrue(assemblies)
        for command in assemblies:
            self.assertEqual(command[command.index("--stop_after") + 1], "scaffolding")
            self.assertEqual(command[command.index("--seed") + 1], "7")
        self.assertIn(">s_subsample_100_1\nACGT\n", combined)

    def test_subsampling_seed_is_used(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            ids = [f"r{i}" for i in range(200)]
            forward, reverse = directory / "R1.fq", directory / "R2.fq"
            write_fastq(forward, ids)
            write_fastq(reverse, ids)
            subsets = {}
            for seed in (1, 1, 2):
                out1, out2 = directory / f"o{seed}_1.fq", directory / f"o{seed}_2.fq"
                MATAM_16S.subsample_pair(str(forward), str(reverse), str(out1), str(out2), 50, seed)
                subsets.setdefault(seed, []).append(read_ids(out1))
            self.assertEqual(subsets[1][0], subsets[1][1])
            self.assertNotEqual(subsets[1][0], subsets[2][0])


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


class MateNamingTests(unittest.TestCase):
    """MarkerMAG matches mates by the naming each sample uses, detected from
    its reads; the read IDs themselves are never changed."""

    NAMINGS = (
        ("shared", ["p1", "p2"], ["p1", "p2"], None),
        ("slash", ["p1/1", "p2/1"], ["p1/2", "p2/2"], "/"),
        ("dot", ["Heavy_1.1", "Heavy_2.1"], ["Heavy_1.2", "Heavy_2.2"], "."),
        ("underscore", ["p1_1", "p2_1"], ["p1_2", "p2_2"], "_"),
        ("colon", ["p1:1", "p2:1"], ["p1:2", "p2:2"], ":"),
        ("sra spots", ["SRR1.1", "SRR1.2"], ["SRR1.1", "SRR1.2"], None),
    )

    def test_detects_each_naming_in_fastq_fasta_and_gzip(self):
        for label, r1, r2, separator in self.NAMINGS:
            for suffix in (".fq.gz", ".fq", ".fasta", ".reads"):
                with self.subTest(naming=label, suffix=suffix), \
                        tempfile.TemporaryDirectory() as directory:
                    directory = Path(directory)
                    paths = []
                    for mate, ids in (("R1", r1), ("R2", r2)):
                        path = directory / f"s_{mate}{suffix}"
                        if suffix in (".fasta",):
                            path.write_text("".join(f">{i}\nACGT\n" for i in ids))
                        else:
                            write_fastq(path, ids)
                        paths.append(path)
                    self.assertEqual(PAIRING.detect_mate_separator(*paths), separator)

    def test_mixed_or_out_of_step_naming_is_an_error(self):
        for r1, r2, message in ((["p1.1", "p2"], ["p1.2", "p2"], "Mixed mate naming"),
                                (["p1.1"], ["p2.2"], "out of step"),
                                (["p1", "p2"], ["p1"], "Cannot detect")):
            with self.subTest(message=message), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                write_fastq(directory / "r1.fq", r1)
                write_fastq(directory / "r2.fq", r2)
                with self.assertRaisesRegex(ValueError, message):
                    PAIRING.detect_mate_separator(directory / "r1.fq", directory / "r2.fq")

    def test_pair_id_strips_only_the_detected_suffix(self):
        cases = (("Heavy_1691.1", ".", "Heavy_1691"), ("Heavy_1691.2 desc", ".", "Heavy_1691"),
                 ("SRR1.1", None, "SRR1.1"), ("p1/2", "/", "p1"), ("p1_1", "/", "p1_1"),
                 ("p1_1", "_", "p1"), ("p1.3", ".", "p1.3"), ("sample.1", "/", "sample.1"))
        for read_id, separator, expected in cases:
            with self.subTest(read_id=read_id, separator=separator):
                self.assertEqual(PAIRING.canonical_read_id(read_id, separator), expected)

    def test_sam_mates_of_dot_named_pair_share_the_pair_id(self):
        # bowtie2 keeps .1/.2 in QNAME (it trims only /1 /2).
        self.assertEqual(PAIRING.sam_read_id_and_mate(["Heavy_1691.1", "99"], "."),
                         ("Heavy_1691", "1"))
        self.assertEqual(PAIRING.sam_read_id_and_mate(["Heavy_1691.2", "147"], "."),
                         ("Heavy_1691", "2"))

    def test_reads_extracted_by_pair_id_keep_their_names(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "r1.fasta"
            source.write_text(">Heavy_1.1\nACGT\n>Heavy_2.1\nACGT\n")
            output = directory / "out.fasta"
            written = PAIRING.extract_reads_by_pair_id(
                source, {"Heavy_2"}, output, sequence_format="fasta", mate_separator=".")
            self.assertEqual(written, 1)
            self.assertEqual(output.read_text(), ">Heavy_2.1\nACGT\n")

    @staticmethod
    def pair_records(name1, name2, pos, seq="A" * 150):
        return [[name1, "99", "r16s", str(pos), "42", "150=", "=", str(pos + 200), "350", seq, "*"],
                [name2, "147", "r16s", str(pos + 200), "42", "150=", "=", str(pos), "-350", seq, "*"]]

    @staticmethod
    def write_sam(path, records):
        Path(path).write_text("@HD\tVN:1.0\tSO:queryname\n"
                              + "".join("\t".join(r) + "\n" for r in records))

    def parse_16s(self, sam, out, separator):
        link_16s = importlib.import_module("MarkerMAG.link_16s")
        link_16s.parse_sam16s_worker([str(sam), str(out), 50, 2, {"r16s": 1500}, separator])
        return Path(out).read_text().splitlines()

    def test_linking_sees_both_mates_of_each_naming(self):
        """link_16s.parse_sam16s_worker on one pair named three ways: the same
        record once the mates are matched by the detected naming. Under the
        old /1-/2-only rule, .1/.2 mates were parsed as separate reads."""
        both_mapped = ["p1", "1", "0", "0", "1", "r16s:::101:150=", "r16s:::301:150=",
                       "na", "na", "r16s:150=,150="]
        for label, names, separator, expected in (
                ("shared", ("p1", "p1"), None, 1), ("slash", ("p1/1", "p1/2"), "/", 1),
                ("dot", ("p1.1", "p1.2"), ".", 1), ("dot, old rule", ("p1.1", "p1.2"), "/", 2)):
            with self.subTest(naming=label), tempfile.TemporaryDirectory() as directory:
                sam = Path(directory) / "reads.sam"
                self.write_sam(sam, self.pair_records(*names, 101))
                lines = self.parse_16s(sam, Path(directory) / "out.mp", separator)
                # The only (so also the last) pair is written: it used to be lost.
                self.assertEqual(len(lines), expected)
                if expected == 1:
                    self.assertEqual(lines[0].split("\t"), both_mapped)

    def test_sam_is_split_between_pairs_only(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            records = []
            for number in range(1, 8):
                records += self.pair_records(f"Heavy_{number}.1", f"Heavy_{number}.2", 101)
                records.append(records[-1][:])  # a second alignment of mate 2
            self.write_sam(directory / "reads.sam", records)
            for lines_per_piece in (1, 2, 3, 4, 100):
                with self.subTest(lines_per_piece=lines_per_piece):
                    folder = directory / f"split_{lines_per_piece}"
                    folder.mkdir()
                    pieces = PAIRING.split_sam_by_pair(directory / "reads.sam", folder,
                                                       lines_per_piece, ".")
                    pairs_by_piece = []
                    body = []
                    for piece in pieces:
                        lines = Path(piece).read_text().splitlines()
                        body += lines
                        pairs_by_piece.append({PAIRING.canonical_read_id(l.split("\t")[0], ".")
                                               for l in lines})
                    self.assertEqual(body, ["\t".join(r) for r in records])  # no header, in order
                    for a, b in zip(pairs_by_piece, pairs_by_piece[1:]):
                        self.assertFalse(a & b)  # no pair in two pieces
                    if lines_per_piece == 3:  # 3 records per pair: one pair per piece
                        self.assertEqual([len(p) for p in pairs_by_piece], [1] * 7)

    def test_linking_records_do_not_depend_on_the_number_of_pieces(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            records = []
            for number in range(1, 12):
                records += self.pair_records(f"Heavy_{number}.1", f"Heavy_{number}.2", 100 + 7 * number)
            self.write_sam(directory / "reads.sam", records)
            whole = self.parse_16s(directory / "reads.sam", directory / "whole.mp", ".")
            self.assertEqual(len(whole), 11)  # every pair, the last one included
            for lines_per_piece in (1, 3, 5, 8):
                with self.subTest(lines_per_piece=lines_per_piece):
                    folder = directory / f"split_{lines_per_piece}"
                    folder.mkdir()
                    merged = []
                    for number, piece in enumerate(PAIRING.split_sam_by_pair(
                            directory / "reads.sam", folder, lines_per_piece, ".")):
                        merged += self.parse_16s(piece, folder / f"{number}.mp", ".")
                    self.assertEqual(merged, whole)


if __name__ == "__main__":
    unittest.main()
