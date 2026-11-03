"""Byte comparisons against executed original Perl, without mocking BEDTools.

The intersection fixtures are explicit BED4+BED5 rows. They test the join
operation independently of interval intersection, which requires real BEDTools.
"""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
CURRENT = Path(os.environ.get("PEREGRINE_TEST_REPO", ROOT.parent)).resolve()
ORIGINAL = ROOT / "original"
IMPLEMENTATION = CURRENT / "chia_pet" / "ChIA_PET.py"
SPEC = importlib.util.spec_from_file_location("chia_pet_rewrite", IMPLEMENTATION)
CHIA = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHIA)

P1 = "HUMAN|HGNC=1|UniProtKB=P00001"
P2 = "HUMAN|HGNC=2|UniProtKB=P00002"
A = ("chr1", 1000, 1600)
B = ("chr1", 3000, 3600)
C = ("chr2", 1000, 1600)


@unittest.skipUnless(shutil.which("perl"), "Original Perl interpreter unavailable")
class OriginalEquivalence(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.panther = self.write("panther", f"{P1}\tG1\n{P2}\tG2\n{P1}\tG3\n")
        self.tissue = self.write("tissues", "59\tHCT-116\n55\tK562\n")
        self.mapping = CHIA.load_panther(self.panther)
        self.tissues = CHIA.load_tissues(self.tissue)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, name, content):
        path = self.directory / name
        path.write_text(content, encoding="utf-8")
        return path

    def perl(self, script, *arguments):
        result = subprocess.run(
            ["perl", str(ORIGINAL / script), *map(str, arguments)],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def assert_bytes(self, python_output, perl_output):
        self.assertEqual(python_output.read_bytes(), perl_output.read_bytes())

    @staticmethod
    def overlap(identifier, anchor, first=A, second=B,
                label="HCT116-POLR2A", pvalue="0.01"):
        def region(value):
            return f"{value[0]}:{value[1]}..{value[2]}"
        return "\t".join(map(str, (
            *anchor, identifier, *anchor,
            f"{region(first)}-{region(second)}", f"{label}_{pvalue}",
        ))) + "\n"

    def compare_links(self, enhancer_rows, gene_rows, expected=None):
        enhancer = self.write("enhancers.intersect", "".join(enhancer_rows))
        genes = self.write("genes.intersect", "".join(gene_rows))
        python_output = self.directory / "python.links"
        perl_output = self.directory / "perl.links"
        CHIA.chiagenesenhancers(
            enhancer, genes, python_output, self.mapping, self.tissues,
        )
        self.perl("chiagenesenhancers.pl", enhancer, genes, perl_output,
                  self.panther, self.tissue)
        self.assert_bytes(python_output, perl_output)
        if expected is not None:
            self.assertEqual(python_output.read_text(), expected)
        return python_output, perl_output

    def test_anchor_conversion_cis_trans_counts_and_pvalues(self):
        raw = self.write("mango", (
            "chr1\t1000\t1600\tchr1\t3000\t3600\t0\t1\n"
            "chr1\t1000\t1600\tchr2\t1000\t1600\t999\t0\n"
            "chrX\t1000\t1600\tchrY\t3000\t3600\t2\t1e-15\n"
        ))
        python_output = self.directory / "python.anchors"
        perl_output = self.directory / "perl.anchors"
        CHIA.celltypev2(raw, "HCT116-POLR2A", python_output)
        self.perl("celltypev2.pl", raw, "HCT116-POLR2A", perl_output)
        self.assert_bytes(python_output, perl_output)
        self.assertEqual(len(python_output.read_text().splitlines()), 6)

    def test_promoters_strands_header_mapping_last_row_and_lexical_sort(self):
        genes = self.write("gene_table", "\t".join(CHIA.GENE_HEADER) + "\n" + (
            "G1\t1000\t1400\t900\t1\t2\n"
            "G2\t2000\t2400\t1900\t-1\t10\n"
            "UNMAPPED\t3000\t3400\t2900\t1\tX\n"
            "G1\t1600\t2000\t1500\t1\t2\n"
        ))
        python_output = self.directory / "python.promoters"
        perl_output = self.directory / "perl.promoters"
        CHIA.chiageneTSS(genes, python_output, self.mapping)
        self.perl("chiageneTSS.pl", genes, self.panther, perl_output)
        self.assert_bytes(python_output, perl_output)
        self.assertEqual(python_output.read_text(),
                         "chr10\t2000\t2600\tG2\nchr2\t1000\t1600\tG1\n")

    def test_cis_partner(self):
        self.compare_links([self.overlap("E2", B)], [self.overlap("G1", A)],
                           f"E2\t{P1}\t59\t1\n")

    def test_reverse_direction(self):
        self.compare_links([self.overlap("E2", A)], [self.overlap("G1", B)],
                           f"E2\t{P1}\t59\t1\n")

    def test_trans_partner(self):
        self.compare_links([self.overlap("E2", C, A, C)],
                           [self.overlap("G1", A, A, C)], f"E2\t{P1}\t59\t1\n")

    def test_same_anchor_does_not_link(self):
        self.compare_links([self.overlap("E2", A)], [self.overlap("G1", A)], "")

    def test_explicit_self_interaction_does_link(self):
        self.compare_links([self.overlap("E2", A, A, A)],
                           [self.overlap("G1", A, A, A)], f"E2\t{P1}\t59\t1\n")

    def test_duplicate_interactions_and_overlaps_deduplicate(self):
        self.compare_links([self.overlap("E2", B)] * 3,
                           [self.overlap("G1", A)] * 2, f"E2\t{P1}\t59\t1\n")

    def test_different_proteins_do_not_cross_match(self):
        self.compare_links([self.overlap("E2", B, label="HCT116-CTCF")],
                           [self.overlap("G1", A)], "")

    def test_same_cell_two_proteins_merge_final_link(self):
        self.compare_links(
            [self.overlap("E2", B, label=label)
             for label in ["HCT116-CTCF", "HCT116-POLR2A"]],
            [self.overlap("G1", A, label=label)
             for label in ["HCT116-CTCF", "HCT116-POLR2A"]],
            f"E2\t{P1}\t59\t1\n",
        )

    def test_probability_values_are_not_link_filters_or_output_fields(self):
        self.compare_links(
            [self.overlap("E2", B, pvalue=p) for p in ["1", "0", "0.9", "1e-40"]],
            [self.overlap("G1", A, pvalue=p) for p in ["1", "0", "0.9", "1e-40"]],
            f"E2\t{P1}\t59\t1\n",
        )

    def test_two_ensg_one_panther_duplicate_removed_by_orderlinks(self):
        python_links, perl_links = self.compare_links(
            [self.overlap("E2", B)],
            [self.overlap(gene, A) for gene in ["G1", "G3"]],
            f"E2\t{P1}\t59\t1\n" * 2,
        )
        python_output = self.directory / "python.final"
        perl_output = self.directory / "perl.final"
        CHIA.orderlinks(python_links, python_output)
        self.perl("orderlinks.pl", perl_links, perl_output)
        self.assert_bytes(python_output, perl_output)
        self.assertEqual(python_output.read_text(), f"E2\t{P1}\t59\t1\n")

    def test_orderlinks_lexical_and_duplicates(self):
        raw = self.write("unordered", (
            f"2\t{P2}\t59\t1\n"
            f"10\t{P2}\t59\t1\n"
            f"10\t{P1}\t9\t1\n"
            f"10\t{P1}\t55\t1\n"
            f"10\t{P1}\t55\t1\n"
        ))
        python_output = self.directory / "python.ordered"
        perl_output = self.directory / "perl.ordered"
        CHIA.orderlinks(raw, python_output)
        self.perl("orderlinks.pl", raw, perl_output)
        self.assert_bytes(python_output, perl_output)
        self.assertEqual(python_output.read_text().splitlines()[0], f"10\t{P1}\t55\t1")

    def test_empty_join_and_empty_orderlinks(self):
        python_links, perl_links = self.compare_links([], [], "")
        python_output = self.directory / "python.empty"
        perl_output = self.directory / "perl.empty"
        CHIA.orderlinks(python_links, python_output)
        self.perl("orderlinks.pl", perl_links, perl_output)
        self.assert_bytes(python_output, perl_output)
        self.assertEqual(python_output.read_bytes(), b"")


if __name__ == "__main__":
    unittest.main(verbosity=2)
