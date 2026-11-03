#!/usr/bin/env python3
"""Check the actual PEREGRINE checkout; synthetic tests never use production data.

Run with the checkout's venv Python for --scope full. A real BEDTools executable
is required for ChIA execution and --self-test. No dependencies are installed.
"""
import argparse
import hashlib
import importlib
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


TAD_STAGES = (
    "TADorder", "bedtoolIntersect", "split", "bTADlinks",
    "concatenate_linksbTAD", "tTADlinks", "linksbTAD_tissueReplace",
    "final_concatenate", "cutdowntad", "orderlinks",
)
CHIA = "chia_pet/ChIA_PET.py"
CUTDOWN = "tad/cutdowntad.py"


def run(command, **kwargs):
    subprocess.run([str(value) for value in command], check=True, **kwargs)


def require_files(paths, allow_empty=False):
    missing = [str(path) for path in paths if not path.is_file()
               or (not allow_empty and path.stat().st_size == 0)]
    if missing:
        raise ValueError("Missing or empty required files:\n  " + "\n  ".join(missing))


def check_source(repo, scope):
    scripts = [CHIA, CUTDOWN] if scope == "chia" else [CUTDOWN, "tad/orderlinks.py"]
    if scope == "full":
        scripts = [CHIA, "eQTL/eqtl.py", "heirarchicalTAD/heirarchicalTAD.py"]
        scripts.extend("tad/" + name + ".py" for name in TAD_STAGES)
    shells = [repo / "steps.sh", repo / "tad/tad.sh"]
    require_files([repo / name for name in scripts] + shells)
    for name in scripts:
        path = repo / name
        compile(path.read_bytes(), str(path), "exec")
    for shell in shells:
        run(["bash", "-n", shell])
    print(f"PASS syntax: {len(scripts)} Python entry points and {len(shells)} shell scripts", flush=True)


def check_environment(scope, self_test):
    if sys.version_info < (3, 9):
        raise ValueError("Python 3.9 or newer is required")
    print(f"Python: {sys.executable} ({sys.version.split()[0]})", flush=True)
    if scope == "full":
        for name in ("numpy", "pandas", "pybedtools", "pysam"):
            module = importlib.import_module(name)
            print(f"Dependency: {name} {getattr(module, '__version__', 'import OK')}", flush=True)
    if scope != "final" or self_test:
        executable = shutil.which("bedtools")
        if not executable:
            raise ValueError("Real BEDTools is not on PATH; install/load it before running this check")
        run([executable, "--version"])
    print("PASS environment", flush=True)


def check_inputs(repo, scope, manifest, genes, output):
    if scope in ("chia", "full"):
        require_files([manifest, genes, repo / "pantherGeneList.txt",
                       repo / "tissuetable_10092018.txt", repo / "CREbedDBenhancers_10092018"])
        if scope == "full":
            require_files([repo / path for path in (
                "eQTL/exons_genes.txt", "heirarchicalTAD/UP000005640_9606.idmapping",
                "tad/ENCFF274VJU.bed", "tad/ENCFF588KUZ.bed", "tad/TSSgenesbed",
            )])
            for pattern in ("eQTL/GTEx_Analysis_v7_eQTL/*variant_gene_pairs.txt.gz",
                            "heirarchicalTAD/data/*.bed"):
                if not any(path.is_file() and path.stat().st_size for path in repo.glob(pattern)):
                    raise ValueError("No nonempty inputs match " + str(repo / pattern))
        run([sys.executable, repo / CHIA, output / "chia_pet", "--repo-root", repo,
             "--manifest", manifest, "--genes", genes, "--check-inputs"])
    else:
        # Zero-row evidence/results may be valid. Missing files are not valid.
        require_files([output / path for path in (
            "eQTL/linksDBnumeqtl", "heirarchicalTAD/PSYCHIClinksDB",
            "tad/linksTADtissues", "chia_pet/linksDBchia",
        )], allow_empty=True)
        with (output / "chia_pet/linksDBchia").open(encoding="utf-8") as source:
            for number, line in enumerate(source, 1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 4 or not all(fields) or fields[3] != "1":
                    raise ValueError(f"linksDBchia:{number}: expected four fields ending in assay 1")
    print("PASS production input validation", flush=True)


def self_test(repo):
    """Use the actual uploaded implementation and real BEDTools in a private fixture."""
    require_files([repo / CHIA, repo / CUTDOWN])
    print("SYNTHETIC SELF-TEST: fixture data only; no production files are changed", flush=True)
    with tempfile.TemporaryDirectory(prefix="peregrine synthetic check ") as temporary:
        root = Path(temporary)
        fixture = root / "fixture references"
        (fixture / "chia_pet").mkdir(parents=True)
        (fixture / "pantherGeneList.txt").write_text("PANTHER1\tENSG1\nPANTHER2\tENSG2\nPANTHER3\tENSG3\n")
        (fixture / "tissuetable_10092018.txt").write_text("55\tK562\n59\tHCT116\n")
        (fixture / "CREbedDBenhancers_10092018").write_text(
            "chr1\t2000\t2100\tE_POS\nchr1\t2090\t2110\tE_HALF\n"
            "chr1\t2091\t2111\tE_BELOW\nchr1\t1100\t1200\tE_SAME\n"
            "chr2\t5000\t5100\tE_NEG\nchr1\t6050\t6150\tE_SELF\n")
        (fixture / "chia_pet/hg19genes.txt").write_text(
            "ENSG1\t1600\t1900\t1600\t1\t1\n"
            "ENSG2\t4000\t4000\t3500\t-1\t2\n"
            "ENSG3\t6600\t7000\t6600\t1\t1\n")
        (fixture / "chia_pet/input.mango").write_text(
            "chr1\t1000\t1600\tchr1\t2000\t2100\t4\t0.01\n" * 2
            + "chr2\t4000\t4600\tchr2\t5000\t5100\t3\t0.02\n"
            + "chr1\t6000\t6600\tchr1\t6000\t6600\t2\t0.03\n")
        manifest = fixture / "chia_pet/samples.tsv"
        manifest.write_text("sample_id\tcell_type\ttarget_protein\tmango_file\n"
                            "test_K562\tK562\tPOLR2A\tinput.mango\n"
                            "test_HCT116\tHCT116\tPOLR2A\tinput.mango\n")
        output = root / "synthetic output"
        command = [sys.executable, repo / CHIA, output / "chia_pet", "--repo-root", fixture,
                   "--manifest", manifest, "--genes", fixture / "chia_pet/hg19genes.txt"]
        run(command + ["--check-inputs"])
        run(command)
        expected = "".join("\t".join((enhancer, gene, tissue, "1")) + "\n"
                           for enhancer, gene in (("E_HALF", "PANTHER1"), ("E_NEG", "PANTHER2"),
                                                  ("E_POS", "PANTHER1"), ("E_SELF", "PANTHER3"))
                           for tissue in ("55", "59"))
        result = (output / "chia_pet/linksDBchia").read_text()
        if result != expected:
            raise ValueError("Synthetic ChIA output differs from the expected eight links")
        for folder in ("eQTL", "heirarchicalTAD", "tad"):
            (output / folder).mkdir()
        (output / "eQTL/linksDBnumeqtl").write_text("enhancer\tgene\ttissue\tnumber_of_eQTL\tassay\nE_EQTL\tPE\t1\t1\t2\n")
        (output / "heirarchicalTAD/PSYCHIClinksDB").write_text("E_HIER\tPH\t49\t0.01\t4\n")
        tad_rows = ["E_POS\tPANTHER1\t65\t3\n", "E_NEG\tPANTHER2\t65\t3\n",
                    "E_UNSUPPORTED\tPN\t65\t3\n", "E_EQTL\tPE\t65\t3\n", "E_HIER\tPH\t65\t3\n"]
        (output / "tad/linksTADtissues").write_text("".join(tad_rows))
        run([sys.executable, repo / CUTDOWN, output / "tad"])
        if (output / "tad/selectTAD").read_text() != "".join(tad_rows[3:]):
            raise ValueError("Synthetic TAD filtering without ChIA differs from expected")
        run([sys.executable, repo / CUTDOWN, output / "tad", "--chia-links", output / "chia_pet/linksDBchia"])
        selected = "".join(tad_rows[:2] + tad_rows[3:])
        if (output / "tad/selectTAD").read_text() != selected:
            raise ValueError("Synthetic TAD filtering with ChIA differs from expected")
        order = repo / "tad/orderlinks.py"
        if order.is_file():
            run([sys.executable, order, output / "tad"])
            if (output / "tad/linksDBtad").read_text() != "".join(sorted(set(selected.splitlines(True)))):
                raise ValueError("Synthetic final TAD ordering differs from expected")
            print("PASS synthetic final TAD sorting", flush=True)
        else:
            print("NOT CHECKED: final TAD sorting; tad/orderlinks.py is absent from this source snapshot", flush=True)
        print("PASS synthetic ChIA: 8 links; TAD filter: 2 rows without ChIA, 4 with ChIA", flush=True)
        print("Synthetic linksDBchia SHA256: " + hashlib.sha256(result.encode()).hexdigest(), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--scope", choices=("full", "chia", "final"), default="full")
    parser.add_argument("--syntax-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--check-inputs", action="store_true")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--genes", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.syntax_only and (args.self_test or args.check_inputs):
        parser.error("--syntax-only cannot be combined with runtime/input tests")
    repo = args.repo.resolve()
    manifest = args.manifest.resolve() if args.manifest else repo / "chia_pet/samples.tsv"
    genes = args.genes.resolve() if args.genes else repo / "chia_pet/hg19genes.txt"
    output = args.output.resolve() if args.output else repo / "results"
    try:
        check_source(repo, args.scope)
        if args.syntax_only:
            print("SYNTAX ONLY: dependencies, input data and biological execution were not tested")
            return 0
        check_environment(args.scope, args.self_test)
        if args.check_inputs:
            check_inputs(repo, args.scope, manifest, genes, output)
        if args.self_test:
            self_test(repo)
        if not args.check_inputs:
            print("Production data readiness was not checked; use --check-inputs before a real run", flush=True)
    except (OSError, ValueError, SyntaxError, ImportError, subprocess.CalledProcessError) as error:
        print("CHECK FAILED: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
