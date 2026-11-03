"""Native evidence-filter tests plus explicitly controlled launcher argument tests.

The launcher recorder tests do not run the biological analysis stages.  The
post-Mango end-to-end test uses real BEDTools and skips when it is unavailable.
"""

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
FILES = Path(os.environ.get("PEREGRINE_TEST_REPO", ROOT.parent)).resolve()
CUTDOWN = FILES / "tad" / "cutdowntad.py"
CHIA = FILES / "chia_pet" / "ChIA_PET.py"


class EvidenceFilterIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="chia integration ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.target = self.root / "results with spaces" / "tad"
        self.target.mkdir(parents=True)
        self.eqtl = self.target.parent / "eQTL" / "linksDBnumeqtl"
        self.hierarchy = self.target.parent / "heirarchicalTAD" / "PSYCHIClinksDB"
        self.chia = self.target.parent / "chia_pet" / "linksDBchia"
        for source in (self.eqtl, self.hierarchy, self.chia):
            source.parent.mkdir(parents=True)
            source.write_text("")
        self.tad = self.target / "linksTADtissues"
        self.tad.write_text("")

    def cli(self, chia=None):
        command = [sys.executable, str(CUTDOWN), str(self.target)]
        if chia is not None:
            command.extend(("--chia-links", str(chia)))
        return subprocess.run(command, cwd=self.root, capture_output=True, text=True)

    def selected(self):
        return (self.target / "selectTAD").read_text()

    def test_support_is_union_of_three_evidence_sources(self):
        self.eqtl.write_text(
            "enhancer\tgene\ttissue\tnumber_of_eQTL\tassay\n"
            "Eeqtl\tP1\t1\t2\t2\n"
        )
        self.hierarchy.write_text("Ehierarchy\tP2\t2\t0.01\t4\n")
        self.chia.write_text("Echia\tP3\t3\t1\n")
        kept = "Eeqtl\tP1\t8\t3\nEhierarchy\tP2\t8\t3\nEchia\tP3\t8\t3\n"
        self.tad.write_text(kept + "Eunsupported\tP4\t8\t3\n")
        result = self.cli(self.chia)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.selected(), kept)
        self.assertIn("or ChIA: 3", result.stdout)

    def test_chia_support_ignores_evidence_tissue_and_preserves_tad_assay(self):
        self.chia.write_text("E1\tP1\t99\t1\n")
        self.tad.write_text("E1\tP1\t7\t3\n")
        result = self.cli(self.chia)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.selected(), "E1\tP1\t7\t3\n")

    def test_omitted_option_does_not_auto_load_stale_chia_file(self):
        self.chia.write_text("E1\tP1\t99\t1\n")
        self.tad.write_text("E1\tP1\t7\t3\n")
        result = self.cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.selected(), "")

    def test_missing_explicit_chia_file_fails_before_overwriting_selection(self):
        selected = self.target / "selectTAD"
        selected.write_text("previous complete result\n")
        result = self.cli(self.root / "missing chia file")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing chia file", result.stderr)
        self.assertEqual(selected.read_text(), "previous complete result\n")

    def test_explicit_chia_wrong_assay_is_rejected(self):
        self.chia.write_text("E1\tP1\t99\t4\n")
        result = self.cli(self.chia)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ChIA assay 1", result.stderr)

    def test_tad_order_and_duplicate_rows_are_preserved(self):
        self.eqtl.write_text("E1\tP1\t1\t1\t2\nE2\tP2\t1\t1\t2\n")
        rows = "E2\tP2\t7\t3\nE1\tP1\t8\t3\nE2\tP2\t7\t3\n"
        self.tad.write_text(rows)
        result = self.cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.selected(), rows)

    def test_existing_four_argument_function_call_remains_valid(self):
        specification = importlib.util.spec_from_file_location("cutdown_contract", CUTDOWN)
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        self.eqtl.write_text("E1\tP1\t1\t1\t2\n")
        self.tad.write_text("E1\tP1\t2\t3\n")
        output = self.target / "old function result"
        module.cutdowntad(self.eqtl, self.hierarchy, self.tad, output)
        self.assertEqual(output.read_text(), "E1\tP1\t2\t3\n")


class ControlledLauncherArgumentTests(unittest.TestCase):
    """Recorder-only unit checks; no fake BEDTools or biological execution."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="chia launcher ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.repo = self.root / "repository with spaces"
        (self.repo / "tad").mkdir(parents=True)
        (self.repo / "venv" / "bin").mkdir(parents=True)
        shutil.copyfile(FILES / "tad" / "tad.sh", self.repo / "tad" / "tad.sh")
        self.log = self.root / "argv.jsonl"
        recorder = self.repo / "venv" / "bin" / "python"
        recorder.write_text(
            f"#!{sys.executable}\n"
            "import json, os, sys\n"
            "with open(os.environ['CHIA_TEST_ARGV_LOG'], 'a') as output:\n"
            "    output.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        )
        recorder.chmod(0o755)
        self.environment = dict(os.environ, CHIA_TEST_ARGV_LOG=str(self.log))
        self.chia = self.root / "chia evidence with spaces.tsv"
        self.chia.write_text("E1\tP1\t1\t1\n")

    def launcher(self, extra=()):
        return subprocess.run(
            ["bash", str(self.repo / "tad" / "tad.sh"), "results with spaces/tad", *extra],
            cwd=self.root, env=self.environment, capture_output=True, text=True,
        )

    def recorded(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()]

    def test_chia_path_is_forwarded_only_to_cutdown_stage(self):
        result = self.launcher(("--chia-links", self.chia.name))
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.recorded()
        self.assertEqual(len(calls), 10)
        target = str(self.root / "results with spaces" / "tad")
        for call in calls:
            if Path(call[0]).name == "cutdowntad.py":
                self.assertEqual(call[1:], [target, "--chia-links", str(self.chia)])
            else:
                self.assertEqual(call[1:], [target])

    def test_existing_output_only_invocation_passes_no_chia_option(self):
        result = self.launcher()
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.recorded()
        self.assertEqual(len(calls), 10)
        self.assertTrue(all(len(call) == 2 for call in calls))

    def test_missing_explicit_evidence_fails_before_any_stage(self):
        result = self.launcher(("--chia-links", "missing evidence.tsv"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ChIA links file not found", result.stderr)
        self.assertFalse(self.log.exists())


class RealBedtoolsWorkflowTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("bedtools"), "real BEDTools is not installed")
    def test_post_mango_cli_runs_real_intersections_and_writes_final_links(self):
        with tempfile.TemporaryDirectory(prefix="chia real bedtools ") as temporary:
            repo = Path(temporary)
            (repo / "chia_pet").mkdir()
            (repo / "pantherGeneList.txt").write_text("PANTHER1\tENSG1\n")
            (repo / "tissuetable_10092018.txt").write_text("55\tK562\n")
            (repo / "CREbedDBenhancers_10092018").write_text("chr1\t2000\t2100\tE1\n")
            (repo / "chia_pet" / "hg19genes.txt").write_text("ENSG1\t1600\t1700\t1500\t1\t1\n")
            (repo / "chia_pet" / "sample.mango").write_text(
                "chr1\t1000\t1600\tchr1\t2000\t2100\t4\t0.01\n"
            )
            (repo / "chia_pet" / "samples.tsv").write_text(
                "sample_id\tcell_type\ttarget_protein\tmango_file\n"
                "sample1\tK562\tPOLR2A\tsample.mango\n"
            )
            output = repo / "output with spaces"
            result = subprocess.run(
                [sys.executable, str(CHIA), str(output), "--repo-root", str(repo)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((output / "linksDBchia").read_text(), "E1\tPANTHER1\t55\t1\n")
            self.assertEqual([path.name for path in output.iterdir()], ["linksDBchia"])


if __name__ == "__main__":
    unittest.main()
