"""Wrapper dispatch/argument tests with recorded calls, not biological validation.

The temporary fake stages below only record argv and return a requested status.
They do not produce or validate production PEREGRINE datasets.
"""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


WRAPPER_SOURCE = Path(os.environ.get(
    "PEREGRINE_TEST_WRAPPER",
    Path(__file__).resolve().parents[1] / "run_peregrine.sh",
)).resolve()


class WrapperArguments(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="peregrine-wrapper-unit-")
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name).resolve()
        self.repo = self.base / "repository with spaces"
        self.caller = self.base / "caller with spaces"
        self.repo.mkdir()
        self.caller.mkdir()
        self.wrapper = self.repo / "run_peregrine.sh"
        shutil.copyfile(WRAPPER_SOURCE, self.wrapper)
        self.log = self.base / "calls.jsonl"
        self.environment = dict(os.environ, WRAPPER_CALL_LOG=str(self.log))

        for relative in ("venv/bin", "chia_pet", "tad"):
            (self.repo / relative).mkdir(parents=True)
        (self.repo / "venv/bin/python").symlink_to(sys.executable)

        recorder = """import json, os, pathlib, sys
program = __PROGRAM_NAME__
with open(os.environ['WRAPPER_CALL_LOG'], 'a', encoding='utf-8') as out:
    out.write(json.dumps({'program': program, 'args': sys.argv[1:], 'cwd': str(pathlib.Path.cwd())}) + '\\n')
if os.environ.get('WRAPPER_FAIL_PROGRAM') == program:
    sys.exit(17)
if program == 'verify':
    scope = sys.argv[sys.argv.index('--scope') + 1]
    output = pathlib.Path(sys.argv[sys.argv.index('--output') + 1])
    if scope == 'final' and not (output / 'chia_pet/linksDBchia').is_file():
        sys.exit(19)
"""
        for relative, program in (
            ("verify_peregrine.py", "verify"),
            ("chia_pet/ChIA_PET.py", "chia"),
            ("tad/cutdowntad.py", "cutdowntad"),
            ("tad/orderlinks.py", "orderlinks"),
            ("record_full.py", "full"),
        ):
            (self.repo / relative).write_text(recorder.replace("__PROGRAM_NAME__", repr(program)), encoding="utf-8")

        (self.repo / "steps.sh").write_text(
            "#!/usr/bin/env bash\nexec "
            + shlex.quote(sys.executable)
            + " "
            + shlex.quote(str(self.repo / "record_full.py"))
            + ' "$@"\n',
            encoding="utf-8",
        )

    def invoke(self, *args, failure=None):
        environment = dict(self.environment)
        if failure:
            environment["WRAPPER_FAIL_PROGRAM"] = failure
        return subprocess.run(
            ["bash", str(self.wrapper), *args], cwd=self.caller,
            env=environment, text=True, capture_output=True, check=False,
        )

    def calls(self):
        if not self.log.exists():
            return []
        return [json.loads(line) for line in self.log.read_text(encoding="utf-8").splitlines()]

    def clear_calls(self):
        if self.log.exists():
            self.log.unlink()

    def make_final_output(self):
        output = self.caller / "new output"
        (output / "chia_pet").mkdir(parents=True)
        # A successful zero-link result is valid; only existence is required here.
        (output / "chia_pet/linksDBchia").touch()
        return output

    def test_full_delegates_existing_steps_with_default_manifest(self):
        result = self.invoke("new output")
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual([x["program"] for x in calls], ["verify", "full"])
        output = str(self.caller / "new output")
        manifest = str(self.repo / "chia_pet/samples.tsv")
        genes = str(self.repo / "chia_pet/hg19genes.txt")
        self.assertEqual(calls[0]["args"], [
            "--repo", str(self.repo), "--scope", "full", "--output", output,
            "--manifest", manifest, "--genes", genes, "--check-inputs",
        ])
        self.assertEqual(calls[1]["args"], [output, "--chia-manifest", manifest, "--chia-genes", genes])
        self.assertTrue(all(x["cwd"] == str(self.repo) for x in calls))

    def test_only_chia_resolves_custom_paths_before_cd_and_uses_fallback(self):
        (self.repo / "venv/bin/python").unlink()
        result = self.invoke(
            "new output", "--only-chia", "--chia-manifest", "input list.tsv",
            "--chia-genes", "gene annotation.tsv",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual([x["program"] for x in calls], ["verify", "chia"])
        self.assertEqual(calls[0]["args"][3], "chia")
        self.assertIn("--check-inputs", calls[0]["args"])
        self.assertEqual(calls[1]["args"], [
            str(self.caller / "new output/chia_pet"),
            "--manifest", str(self.caller / "input list.tsv"),
            "--genes", str(self.caller / "gene annotation.tsv"),
        ])
        self.assertEqual(result.stdout, "ChIA-PET done\n")

    def test_final_uses_existing_chia_output_without_manifest_arguments(self):
        output = self.make_final_output()
        (self.repo / "venv/bin/python").unlink()
        result = self.invoke(str(output), "--final-only")
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual([x["program"] for x in calls], ["verify", "cutdowntad", "orderlinks"])
        self.assertEqual(calls[0]["args"], [
            "--repo", str(self.repo), "--scope", "final", "--output", str(output), "--check-inputs",
        ])
        self.assertEqual(calls[1]["args"], [
            str(output / "tad"), "--chia-links", str(output / "chia_pet/linksDBchia"),
        ])
        self.assertEqual(calls[2]["args"], [str(output / "tad")])
        self.assertEqual(result.stdout, "cutdowntad done\norderlinks done\nTAD done\n")

    def test_check_inputs_never_runs_biological_stages(self):
        output = self.make_final_output()
        for mode, scope in ((None, "full"), ("--only-chia", "chia"), ("--final-only", "final")):
            with self.subTest(mode=mode):
                self.clear_calls()
                args = [str(output), "--check-inputs"]
                if mode:
                    args.append(mode)
                result = self.invoke(*args)
                self.assertEqual(result.returncode, 0, result.stderr)
                calls = self.calls()
                self.assertEqual(len(calls), 1)
                self.assertEqual(calls[0]["program"], "verify")
                self.assertEqual(calls[0]["args"][3], scope)
                self.assertIn("--check-inputs", calls[0]["args"])

    def test_invalid_options_fail_before_any_child_runs(self):
        cases = (
            (), ("",), ("--only-chia",), ("out", "--unknown"),
            ("out", "--only-chia", "--final-only"),
            ("out", "--only-chia", "--only-chia"),
            ("out", "--final-only", "--final-only"),
            ("out", "--chia-manifest"), ("out", "--chia-genes", ""),
            ("out", "--chia-manifest", "--only-chia"),
            ("out", "--chia-manifest", "a", "--chia-manifest", "b"),
            ("out", "--chia-genes", "a", "--chia-genes", "b"),
            ("out", "--check-inputs", "--check-inputs"),
            ("out", "--final-only", "--chia-manifest", "a"),
            ("out", "--final-only", "--chia-genes", "a"),
        )
        for args in cases:
            with self.subTest(args=args):
                self.clear_calls()
                result = self.invoke(*args)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.calls(), [])

    def test_full_requires_the_existing_environment(self):
        (self.repo / "venv/bin/python").unlink()
        result = self.invoke("out")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Full runs require the existing environment", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_failed_verifier_prevents_all_stage_execution(self):
        for mode in (None, "--only-chia", "--final-only"):
            with self.subTest(mode=mode):
                self.clear_calls()
                args = ["out"] + ([mode] if mode else [])
                result = self.invoke(*args, failure="verify")
                self.assertEqual(result.returncode, 17)
                self.assertEqual([x["program"] for x in self.calls()], ["verify"])
                self.assertNotIn("done", result.stdout)

    def test_final_requires_verifier_approval_of_existing_chia_file(self):
        result = self.invoke("missing output", "--final-only")
        self.assertEqual(result.returncode, 19)
        self.assertEqual([x["program"] for x in self.calls()], ["verify"])

    def test_child_failures_propagate_and_stop_the_sequence(self):
        output = self.make_final_output()
        for mode, failure, expected_calls, stdout in (
            (None, "full", ["verify", "full"], ""),
            ("--only-chia", "chia", ["verify", "chia"], ""),
            ("--final-only", "cutdowntad", ["verify", "cutdowntad"], ""),
            ("--final-only", "orderlinks", ["verify", "cutdowntad", "orderlinks"], "cutdowntad done\n"),
        ):
            with self.subTest(mode=mode, failure=failure):
                self.clear_calls()
                args = [str(output)] + ([mode] if mode else [])
                result = self.invoke(*args, failure=failure)
                self.assertEqual(result.returncode, 17)
                self.assertEqual([x["program"] for x in self.calls()], expected_calls)
                self.assertEqual(result.stdout, stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
