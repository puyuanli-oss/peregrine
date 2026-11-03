"""CURRENT-upload checks against original Perl. All inputs are synthetic tests.

Real intersection tests require an actual BEDTools executable and never use a
stub. Original Perl scripts are extracted from the user's original source zip.
"""

from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent
CURRENT = Path(os.environ.get("PEREGRINE_TEST_REPO", ROOT.parent)).resolve()
ORIGINAL = ROOT / "original"
P1 = "HUMAN|HGNC=1|UniProtKB=P00001"
P2 = "HUMAN|HGNC=2|UniProtKB=P00002"


@unittest.skipUnless(shutil.which("perl"), "requires original Perl")
class OriginalTadFilterEquivalence(unittest.TestCase):
    def test_enabled_and_disabled_chia_match_original_perl_bytes(self):
        with tempfile.TemporaryDirectory(prefix="synthetic tad equivalence ") as tmp:
            root = Path(tmp)
            tad = root / "tad"
            tad.mkdir()
            eqtl = root / "eQTL" / "linksDBnumeqtl"
            hierarchy = root / "heirarchicalTAD" / "PSYCHIClinksDB"
            chia = root / "chia_pet" / "linksDBchia"
            for path in (eqtl, hierarchy, chia):
                path.parent.mkdir()
            eqtl.write_text("enhancer\tgene\ttissue\tnumber_of_eQTL\tassay\n"
                            f"Eeqtl\t{P1}\t1\t2\t2\n")
            hierarchy.write_text(f"Ehierarchy\t{P2}\t2\t0.01\t4\n")
            chia.write_text(f"Echia\t{P1}\t99\t1\n")
            candidates = tad / "linksTADtissues"
            expected_without = (f"Eeqtl\t{P1}\t8\t3\n"
                                f"Ehierarchy\t{P2}\t7\t3\n"
                                f"Eeqtl\t{P1}\t8\t3\n")
            chia_row = f"Echia\t{P1}\t6\t3\n"
            candidates.write_text(expected_without + chia_row + f"UNSUPPORTED\t{P2}\t8\t3\n")
            empty = root / "empty.chia"
            empty.write_text("")
            for enabled in (False, True):
                with self.subTest(chia_enabled=enabled):
                    perl_result = root / f"perl-{enabled}.tsv"
                    perl = subprocess.run([
                        "perl", str(ORIGINAL / "cutdowntad.pl"),
                        str(chia if enabled else empty), str(eqtl), str(hierarchy),
                        str(candidates), str(perl_result)], capture_output=True, text=True)
                    self.assertEqual(perl.returncode, 0, perl.stderr)
                    command = [sys.executable, str(CURRENT / "tad/cutdowntad.py"), str(tad)]
                    if enabled:
                        command.extend(("--chia-links", str(chia)))
                    python = subprocess.run(command, capture_output=True, text=True)
                    self.assertEqual(python.returncode, 0, python.stderr)
                    actual = (tad / "selectTAD").read_bytes()
                    self.assertEqual(actual, perl_result.read_bytes())
                    self.assertEqual(actual.decode(), expected_without + (chia_row if enabled else ""))


@unittest.skipUnless(shutil.which("bedtools") and shutil.which("perl"),
                     "requires real BEDTools and original Perl")
class RealBedtoolsOriginalWorkflowEquivalence(unittest.TestCase):
    def test_every_post_mango_stage_matches_original_perl_bytes(self):
        with tempfile.TemporaryDirectory(prefix="synthetic original ChIA e2e ") as tmp:
            repo = Path(tmp)
            chia_dir = repo / "chia_pet"
            chia_dir.mkdir()
            panther = repo / "pantherGeneList.txt"
            panther.write_text(f"{P1}\tG1\n{P2}\tG2\n{P1}\tG3\n")
            tissues = repo / "tissuetable_10092018.txt"
            tissues.write_text("59\tHCT-116\n55\tK562\n")
            enhancers = repo / "CREbedDBenhancers_10092018"
            enhancers.write_text("chr1\t3000\t3600\t10\nchr2\t5000\t5100\t2\n"
                                 "chr1\t1000\t1600\tSELF\nchr1\t800\t2800\tWIDE\n")
            genes = chia_dir / "hg19genes.txt"
            genes.write_text("G1\t1600\t2000\t1400\t1\t1\n"
                             "G2\t4000\t4300\t3800\t-1\t1\n"
                             "G3\t1600\t2000\t1400\t1\t1\n")
            mango = chia_dir / "synthetic interactions.mango"
            mango.write_text("chr1\t1100\t1200\tchr1\t3050\t3150\t0\t0.9\n"
                             "chr1\t4100\t4200\tchr2\t5020\t5060\t999\t1\n"
                             "chr1\t1100\t1200\tchr1\t1100\t1200\t1\t0\n")
            samples = (("sampleA", "HCT-116", "POLR2A"), ("sampleB", "K562", "CTCF"))
            manifest = chia_dir / "samples.tsv"
            manifest.write_text("sample_id\tcell_type\ttarget_protein\tmango_file\n" + "".join(
                f"{sample}\t{cell}\t{protein}\t{mango.name}\n"
                for sample, cell, protein in samples))
            output = repo / "python results"
            completed = subprocess.run([
                sys.executable, str(CURRENT / "chia_pet/ChIA_PET.py"), str(output),
                "--repo-root", str(repo), "--keep-intermediates"],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            perl_dir = repo / "original results"
            perl_dir.mkdir()

            def perl(script, *arguments):
                result = subprocess.run(["perl", str(ORIGINAL / script), *map(str, arguments)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)

            def intersect(features, anchors, destination):
                with destination.open("wb") as out:
                    result = subprocess.run([
                        shutil.which("bedtools"), "intersect", "-wa", "-wb", "-F", "0.5",
                        "-f", "0.5", "-e", "-a", str(features), "-b", str(anchors)],
                        stdout=out, stderr=subprocess.PIPE)
                self.assertEqual(result.returncode, 0, result.stderr.decode())

            promoters = perl_dir / "TSSgenesbed"
            perl("chiageneTSS.pl", genes, panther, promoters)
            self.assertEqual(promoters.read_bytes(),
                             (output / "intermediates/TSSgenesbed").read_bytes())
            combined = b""
            for sample, cell, protein in samples:
                anchors = perl_dir / f"{sample}.anchors"
                enhancer_overlap = perl_dir / f"{sample}.enhancers"
                gene_overlap = perl_dir / f"{sample}.genes"
                links = perl_dir / f"{sample}.links"
                label = f"{cell.replace('-', '')}-{protein}"
                perl("celltypev2.pl", mango, label, anchors)
                intersect(enhancers, anchors, enhancer_overlap)
                intersect(promoters, anchors, gene_overlap)
                perl("chiagenesenhancers.pl", enhancer_overlap, gene_overlap, links,
                     panther, tissues)
                for expected_path, name in ((anchors, "anchors.bed"),
                                            (enhancer_overlap, "enhancers.intersect"),
                                            (gene_overlap, "genes.intersect"), (links, "links")):
                    self.assertEqual(expected_path.read_bytes(),
                                     (output / "intermediates" / sample / name).read_bytes(), name)
                combined += links.read_bytes()
            raw = perl_dir / "linkschia"
            raw.write_bytes(combined)
            final = perl_dir / "linksDBchia"
            perl("orderlinks.pl", raw, final)
            self.assertEqual(final.read_bytes(), (output / "linksDBchia").read_bytes())
            rows = final.read_text().splitlines()
            self.assertEqual(len(rows), 8)
            self.assertIn(f"WIDE\t{P1}\t59\t1", rows)
            self.assertIn(f"2\t{P2}\t55\t1", rows)
            self.assertIn("ChIA-PET linksDBchia done: 8 links", completed.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
