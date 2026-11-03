"""Verifier preflight/failure tests. Synthetic fixtures; no substitute BEDTools."""

import contextlib
import importlib.util
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parent
CURRENT = Path(os.environ.get("PEREGRINE_TEST_REPO", ROOT.parent)).resolve()
VERIFIER = Path(os.environ.get("PEREGRINE_TEST_VERIFIER", CURRENT / "verify_peregrine.py"))
SPEC = importlib.util.spec_from_file_location("peregrine_verifier_under_test", VERIFIER)
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


class VerifierTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="synthetic verifier check ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.repo = self.root / "checkout with spaces"
        for folder in ("chia_pet", "tad"):
            (self.repo / folder).mkdir(parents=True)
        for path in ("chia_pet/ChIA_PET.py", "tad/cutdowntad.py", "steps.sh", "tad/tad.sh"):
            shutil.copyfile(CURRENT / path, self.repo / path)
        (self.repo / "pantherGeneList.txt").write_text("P1\tENSG1\n")
        (self.repo / "tissuetable_10092018.txt").write_text("55\tK562\n")
        (self.repo / "CREbedDBenhancers_10092018").write_text("chr1\t3000\t3100\tE1\n")
        self.genes = self.repo / "chia_pet/hg19genes.txt"
        self.genes.write_text("ENSG1\t1600\t1800\t1600\t1\t1\n")
        (self.repo / "chia_pet/input.mango").write_text(
            "chr1\t1000\t1600\tchr1\t3000\t3100\t3\t0.01\n")
        self.manifest = self.repo / "chia_pet/samples.tsv"
        self.manifest.write_text("sample_id\tcell_type\ttarget_protein\tmango_file\n"
                                 "SYNTHETIC\tK562\tPOLR2A\tinput.mango\n")
        self.output = self.root / "output must not be created by preflight"

    def cli(self, *extra):
        return subprocess.run([sys.executable, str(VERIFIER), "--repo", str(self.repo),
                               "--scope", "chia", *extra], capture_output=True, text=True)

    def test_syntax_only_needs_no_biological_execution_or_output(self):
        result = self.cli("--syntax-only")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PASS syntax: 2 Python entry points", result.stdout)
        self.assertIn("SYNTAX ONLY", result.stdout)
        self.assertFalse(self.output.exists())

    def test_invalid_python_reports_clean_preflight_failure(self):
        (self.repo / "chia_pet/ChIA_PET.py").write_text("def invalid(:\n")
        result = self.cli("--syntax-only")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CHECK FAILED:", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_invalid_shell_is_detected(self):
        (self.repo / "steps.sh").write_text("#!/bin/bash\nif then\n")
        result = self.cli("--syntax-only")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CHECK FAILED:", result.stderr)

    def test_missing_source_file_is_detected(self):
        (self.repo / "tad/cutdowntad.py").unlink()
        result = self.cli("--syntax-only")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cutdowntad.py", result.stderr)

    def test_full_source_check_requires_all_real_entry_points(self):
        with self.assertRaisesRegex(ValueError, "eQTL/eqtl.py"):
            VERIFY.check_source(self.repo, "full")

    def test_source_only_and_runtime_flags_cannot_be_combined(self):
        result = self.cli("--syntax-only", "--check-inputs")
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot be combined", result.stderr)

    def test_missing_bedtools_is_reported_without_any_substitute(self):
        with mock.patch.object(VERIFY.shutil, "which", return_value=None):
            with self.assertRaisesRegex(ValueError, "Real BEDTools is not on PATH"):
                VERIFY.check_environment("chia", False)

    def test_real_current_chia_input_preflight_creates_no_output(self):
        VERIFY.check_inputs(self.repo, "chia", self.manifest, self.genes, self.output)
        self.assertFalse(self.output.exists())

    def test_missing_mango_file_propagates_real_input_validation_failure(self):
        (self.repo / "chia_pet/input.mango").unlink()
        with self.assertRaises(subprocess.CalledProcessError):
            VERIFY.check_inputs(self.repo, "chia", self.manifest, self.genes, self.output)
        self.assertFalse(self.output.exists())

    def test_final_inputs_allow_zero_rows_but_reject_wrong_chia_assay(self):
        required = ("eQTL/linksDBnumeqtl", "heirarchicalTAD/PSYCHIClinksDB",
                    "tad/linksTADtissues", "chia_pet/linksDBchia")
        for name in required:
            path = self.output / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("")
        VERIFY.check_inputs(self.repo, "final", self.manifest, self.genes, self.output)
        (self.output / "chia_pet/linksDBchia").write_text("E1\tP1\t55\t4\n")
        with self.assertRaisesRegex(ValueError, "assay 1"):
            VERIFY.check_inputs(self.repo, "final", self.manifest, self.genes, self.output)

    def test_full_input_preflight_rejects_missing_or_empty_raw_inputs(self):
        for name in ("eQTL/exons_genes.txt", "heirarchicalTAD/UP000005640_9606.idmapping",
                     "tad/ENCFF274VJU.bed", "tad/ENCFF588KUZ.bed", "tad/TSSgenesbed"):
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("synthetic file-presence fixture only\n")
        with self.assertRaisesRegex(ValueError, "No nonempty inputs match"):
            VERIFY.check_inputs(self.repo, "full", self.manifest, self.genes, self.output)
        self.assertFalse(self.output.exists())

    def test_nonzero_real_subprocess_is_not_silently_accepted(self):
        with self.assertRaises(subprocess.CalledProcessError) as caught:
            VERIFY.run([sys.executable, "-c", "raise SystemExit(7)"])
        self.assertEqual(caught.exception.returncode, 7)


if __name__ == "__main__":
    unittest.main(verbosity=2)
