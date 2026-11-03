"""Native CLI validation; no BEDTools substitute or mocked subprocesses."""

import ast
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
CURRENT = Path(os.environ.get("PEREGRINE_TEST_REPO", ROOT.parent)).resolve()
SOURCE = CURRENT / "chia_pet" / "ChIA_PET.py"


class CliValidationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="chia cli validation ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.repo = self.root / "repository with spaces"
        self.repo.mkdir()
        self.cwd = self.root / "different working directory"
        self.cwd.mkdir()
        self.output = self.cwd / "output with spaces"
        (self.repo / "pantherGeneList.txt").write_text("P1\tENSG1\n")
        (self.repo / "tissuetable_10092018.txt").write_text("55\tK562\n")
        (self.repo / "CREbedDBenhancers_10092018").write_text(
            "chr1\t3000\t3100\tE1\n"
        )
        self.genes = self.repo / "genes with spaces.tsv"
        self.genes.write_text("ENSG1\t1000\t1500\t1000\t1\t1\n")
        self.mango = self.repo / "sample with spaces.mango"
        self.mango.write_text(
            "chr1\t500\t1000\tchr1\t3000\t3100\t1\t0.1\n"
        )
        self.manifest = self.repo / "samples with spaces.tsv"
        self.header = "sample_id\tcell_type\ttarget_protein\tmango_file\n"
        self.sample = "sample1\tK562\tPOLR2A\tsample with spaces.mango\n"
        self.manifest.write_text(self.header + self.sample)

    def cli(self, *, check=True, manifest=None, env=None):
        command = [
            sys.executable,
            str(SOURCE),
            "output with spaces",
            "--repo-root",
            str(self.repo),
            "--manifest",
            str(self.manifest if manifest is None else manifest),
            "--genes",
            str(self.genes),
        ]
        if check:
            command.append("--check-inputs")
        return subprocess.run(
            command,
            cwd=self.cwd,
            env=env,
            capture_output=True,
            text=True,
        )

    def test_preflight_other_cwd_and_paths_with_spaces(self):
        result = self.cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ChIA-PET inputs OK: 1 sample(s)", result.stdout)
        self.assertFalse(self.output.exists())

    def test_missing_manifest_fails_without_output(self):
        missing = self.repo / "missing manifest.tsv"
        result = self.cli(manifest=missing)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ERROR:", result.stderr)
        self.assertIn("missing manifest.tsv", result.stderr)
        self.assertFalse(self.output.exists())

    def test_duplicate_sample_id_is_rejected(self):
        self.manifest.write_text(self.header + self.sample + self.sample)
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("duplicate sample_id", result.stderr)
        self.assertFalse(self.output.exists())

    def test_malformed_mango_is_rejected(self):
        self.mango.write_text("chr1\t500\t1000\tchr1\t3000\t3100\t1\n")
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("expected 8 tab-separated fields", result.stderr)
        self.assertFalse(self.output.exists())

    def test_missing_bedtools_exits_nonzero_without_output(self):
        empty_bin = self.root / "empty bin"
        empty_bin.mkdir()
        environment = dict(os.environ, PATH=str(empty_bin))
        result = self.cli(check=False, env=environment)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Required command not found: bedtools", result.stderr)
        self.assertFalse(self.output.exists())

    def test_python39_grammar_and_import_without_input_io(self):
        ast.parse(SOURCE.read_text(), feature_version=(3, 9))
        code = (
            "import importlib.util; "
            f"spec=importlib.util.spec_from_file_location('chia_import_check', {str(SOURCE)!r}); "
            "module=importlib.util.module_from_spec(spec); "
            "spec.loader.exec_module(module); "
            "print('IMPORT OK')"
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=self.cwd,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "IMPORT OK")
        self.assertEqual(list(self.cwd.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
