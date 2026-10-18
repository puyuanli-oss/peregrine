import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from eQTL import eqtl


class TestEqtl(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.input = self.root / "input.tsv"
        self.output = self.root / "output.tsv"

    def test_trim_preserves_tissue_identity(self):
        self.input.write_text(
            "variant_id\tgene_id\tpval_nominal\n"
            "10_100_A_G_b37\tENSG000001.2\t0.05\n"
        )
        for tissue in ("Liver", "Lung"):
            with self.subTest(tissue=tissue):
                eqtl.trim(self.input, self.output, tissue)
                with self.output.open() as stream:
                    records = list(csv.DictReader(stream, delimiter="\t"))
                self.assertEqual(records[0]["tissue"], tissue)

    def test_trim_rejects_missing_columns(self):
        self.input.write_text("variant_id\n10_100_A_G_b37\n")
        with self.assertRaisesRegex(ValueError, "gene_id"):
            eqtl.trim(self.input, self.output, "Liver")

    def test_full_chromosome_names(self):
        chromosomes = ("1", "10", "22", "X", "Y")
        self.input.write_text(
            "variant\tgene\tpval_nominal\ttissue\n"
            + "".join(
                f"{chrom}_100_A_G_b37\tENSG000001.2\t0.05\tLiver\n"
                for chrom in chromosomes
            )
        )
        eqtl.eqtl_process(self.input, self.output)
        fields = [
            line.split("\t")[:3]
            for line in self.output.read_text().splitlines()
        ]
        self.assertEqual(
            fields,
            [[f"chr{chrom}", "100", "100"] for chrom in chromosomes],
        )

    def test_missing_input_files_fail(self):
        (self.root / "pantherGeneList.txt").write_text(
            "PANTHER1\tENSG000001\n"
        )
        (self.root / "tissuetable_10092018.txt").write_text("T1\tLiver\n")
        result = subprocess.run(
            [
                sys.executable,
                str(Path(eqtl.__file__).resolve()),
                str(self.root / "results"),
            ],
            cwd=self.root,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No eQTL input files", result.stderr)


if __name__ == "__main__":
    unittest.main()