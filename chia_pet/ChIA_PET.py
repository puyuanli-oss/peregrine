#!/usr/bin/env python3
"""Post-Mango ChIA-PET linking, preserving the original PEREGRINE Perl logic."""

import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile


GENE_HEADER = (
    "Gene stable ID",
    "Transcription start site (TSS)",
    "Gene end (bp)",
    "Gene start (bp)",
    "Strand",
    "Chromosome/scaffold name",
)
SAMPLE_HEADER = ("sample_id", "cell_type", "target_protein", "mango_file")


def rows(path, columns, extra=False):
    """Read headerless TSV records; tolerate CRLF, not malformed records."""
    with Path(path).open(encoding="utf-8") as source:
        for number, line in enumerate(source, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < columns or (not extra and len(fields) != columns):
                raise ValueError(f"{path}:{number}: expected {columns} tab-separated fields")
            fields = fields[:columns]
            if not all(fields):
                raise ValueError(f"{path}:{number}: empty required field")
            yield number, fields


def load_panther(path):
    mapping = {}
    for _, fields in rows(path, 2, extra=True):
        mapping[fields[1]] = fields[0]
    if not mapping:
        raise ValueError(f"{path}: no PANTHER mappings")
    return mapping


def load_tissues(path):
    tissues = {}
    for _, fields in rows(path, 2, extra=True):
        tissues[fields[1].replace("-", "")] = fields[0]
    if not tissues:
        raise ValueError(f"{path}: no tissue mappings")
    return tissues


def interval(chromosome, start, end, context):
    try:
        left, right = int(start), int(end)
    except ValueError as error:
        raise ValueError(f"{context}: coordinates must be integers") from error
    if left < 0 or right <= left or any(c.isspace() for c in chromosome):
        raise ValueError(f"{context}: invalid BED interval {chromosome}:{start}-{end}")
    if "-" in chromosome:
        raise ValueError(f"{context}: '-' in chromosome is ambiguous in original anchor names")
    return left, right


def mango_rows(path):
    for number, fields in rows(path, 8):
        context = f"{path}:{number}"
        interval(*fields[:3], context)
        interval(*fields[3:6], context)
        try:
            count, probability = int(fields[6]), float(fields[7])
        except ValueError as error:
            raise ValueError(f"{context}: invalid PET count or p-value") from error
        if count < 0 or not math.isfinite(probability) or not 0 <= probability <= 1:
            raise ValueError(f"{context}: invalid PET count or p-value")
        yield fields


def celltypev2(input_file, cell_type, output_file):
    """Original celltypev2.pl: two BED5 anchors per eight-column interaction."""
    with Path(output_file).open("w", encoding="utf-8") as out:
        for fields in mango_rows(input_file):
            chrom1, start1, end1, chrom2, start2, end2, _, pvalue = fields
            name = f"{chrom1}:{start1}..{end1}-{chrom2}:{start2}..{end2}"
            label = f"{cell_type}_{pvalue}"
            out.write(f"{chrom1}\t{start1}\t{end1}\t{name}\t{label}\n")
            out.write(f"{chrom2}\t{start2}\t{end2}\t{name}\t{label}\n")


def promoter_lines(input_file, panther_mapping):
    genes = {}
    for number, fields in rows(input_file, 6, extra=True):
        if number == 1 and tuple(fields) == GENE_HEADER:
            continue
        gene, tss, _, _, strand, chromosome = fields
        if strand not in ("1", "-1"):
            raise ValueError(f"{input_file}:{number}: strand must be 1 or -1")
        if chromosome.startswith("chr"):
            raise ValueError(f"{input_file}:{number}: gene-table chromosome must omit 'chr'")
        try:
            tss = int(tss)
        except ValueError as error:
            raise ValueError(f"{input_file}:{number}: TSS must be an integer") from error
        start, end = (tss - 600, tss) if strand == "1" else (tss, tss + 600)
        # No clipping, coordinate conversion, or transcript selection is added.
        # As in Perl, the last eligible gene-table row wins.
        genes[gene] = f"chr{chromosome}\t{start}\t{end}\t{gene}"
    result = sorted({genes[gene] for gene in panther_mapping if gene in genes})
    for record in result:
        chrom, start, end, _ = record.split("\t")
        interval(chrom, start, end, str(input_file))
    return result


def chiageneTSS(input_file, output_file, panther_mapping):
    """Original 600-bp, strand-specific, last-record-per-gene promoters."""
    with Path(output_file).open("w", encoding="utf-8") as out:
        for line in promoter_lines(input_file, panther_mapping):
            out.write(line + "\n")


def bedtoolIntersect(gene_file, input_file, output_file):
    """At least half of EITHER the feature or anchor must overlap."""
    command = [
        "bedtools", "intersect", "-wa", "-wb", "-F", "0.5",
        "-f", "0.5", "-e", "-a", str(gene_file), "-b", str(input_file),
    ]
    with Path(output_file).open("w", encoding="utf-8") as out:
        subprocess.run(command, stdout=out, check=True)


def overlap_records(path):
    for number, fields in rows(path, 9):
        _, _, _, identifier, chromosome, start, end, interaction, metadata = fields
        parts = interaction.split("-")
        if len(parts) != 2 or "_" not in metadata:
            raise ValueError(f"{path}:{number}: malformed anchor interaction or cell label")
        label, _ = metadata.rsplit("_", 1)
        if "-" not in label:
            raise ValueError(f"{path}:{number}: expected CELL-PROTEIN label")
        region = f"{chromosome}:{start}..{end}"
        if region not in parts:
            raise ValueError(f"{path}:{number}: overlap anchor is absent from interaction")
        yield identifier, (region, label), (parts[0], label), (parts[1], label)


def chiagenesenhancers(input_file, gene_file, output_file, panther_mapping, tissues):
    """Join promoter-overlapped anchors to enhancer-overlapped partner anchors."""
    interactions, enhancers, genes = {}, {}, {}
    for identifier, region, first, second in overlap_records(input_file):
        interactions.setdefault(first, set()).add(second)
        interactions.setdefault(second, set()).add(first)
        enhancers.setdefault(region, set()).add(identifier)
    for identifier, region, first, second in overlap_records(gene_file):
        interactions.setdefault(first, set()).add(second)
        interactions.setdefault(second, set()).add(first)
        genes.setdefault(identifier, set()).add(region)

    final = set()
    for gene, regions in genes.items():
        for region in regions:
            for partner in interactions.get(region, ()):
                for enhancer in enhancers.get(partner, ()):
                    cell = region[1].split("-", 1)[0]
                    if gene not in panther_mapping or cell not in tissues:
                        raise ValueError(f"Missing PANTHER or tissue mapping for {gene!r}, {cell!r}")
                    final.add((enhancer, gene, cell))
    with Path(output_file).open("w", encoding="utf-8") as out:
        for enhancer, gene, cell in sorted(final):
            out.write(f"{enhancer}\t{panther_mapping[gene]}\t{tissues[cell]}\t1\n")


def orderlinks(input_file, output_file):
    """Deduplicate and sort like orderlinks.pl without an all-links RAM set."""
    output_file = Path(output_file)
    with tempfile.TemporaryDirectory(prefix=".chia-sort-", dir=output_file.parent) as temporary:
        with sqlite3.connect(str(Path(temporary) / "links.sqlite")) as database:
            database.execute(
                "CREATE TABLE links (enhancer TEXT, gene TEXT, tissue TEXT, "
                "PRIMARY KEY (enhancer, gene, tissue)) WITHOUT ROWID"
            )
            def records():
                for number, fields in rows(input_file, 4):
                    if fields[3] != "1":
                        raise ValueError(f"{input_file}:{number}: expected ChIA assay 1")
                    yield tuple(fields[:3])
            database.executemany("INSERT OR IGNORE INTO links VALUES (?, ?, ?)", records())
            with output_file.open("w", encoding="utf-8") as out:
                for enhancer, gene, tissue in database.execute(
                    "SELECT enhancer, gene, tissue FROM links ORDER BY enhancer, gene, tissue"
                ):
                    out.write(f"{enhancer}\t{gene}\t{tissue}\t1\n")


def load_samples(manifest, tissues):
    manifest = Path(manifest)
    samples, seen = [], set()
    with manifest.open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SAMPLE_HEADER:
            raise ValueError(f"{manifest}: expected header " + "\t".join(SAMPLE_HEADER))
        for number, row in enumerate(reader, 2):
            if None in row or not all(row.get(key) for key in SAMPLE_HEADER):
                raise ValueError(f"{manifest}:{number}: expected four nonempty fields")
            sample = row["sample_id"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) or sample in seen:
                raise ValueError(f"{manifest}:{number}: invalid or duplicate sample_id {sample!r}")
            seen.add(sample)
            cell = row["cell_type"].replace("-", "")
            protein = row["target_protein"]
            if not re.fullmatch(r"[A-Za-z0-9.]+", cell) or any(c.isspace() for c in protein):
                raise ValueError(f"{manifest}:{number}: invalid cell_type or target_protein")
            if cell not in tissues:
                raise ValueError(f"{manifest}:{number}: no tissue code for {cell!r}")
            mango = Path(row["mango_file"])
            if not mango.is_absolute():
                mango = manifest.parent / mango
            if not mango.is_file():
                raise FileNotFoundError(f"{manifest}:{number}: Mango input not found: {mango}")
            samples.append((sample, f"{cell}-{protein}", mango))
    if not samples:
        raise ValueError(f"{manifest}: no ChIA-PET samples")
    return samples


def validate_inputs(repo_root, manifest, genes):
    panther = load_panther(repo_root / "pantherGeneList.txt")
    tissues = load_tissues(repo_root / "tissuetable_10092018.txt")
    samples = load_samples(manifest, tissues)
    promoter_lines(genes, panther)
    for number, fields in rows(repo_root / "CREbedDBenhancers_10092018", 4):
        interval(*fields[:3], f"enhancers:{number}")
    for _, _, mango in samples:
        for _ in mango_rows(mango):
            pass
    return panther, tissues, samples


def run_workflow(target_path, repo_root, manifest, genes, keep_intermediates=False):
    panther, tissues, samples = validate_inputs(repo_root, manifest, genes)
    if shutil.which("bedtools") is None:
        raise FileNotFoundError("Required command not found: bedtools")
    target_path.mkdir(parents=True, exist_ok=True)
    temporary = None
    if keep_intermediates:
        work = target_path / "intermediates"
        work.mkdir(parents=True, exist_ok=True)
    else:
        temporary = tempfile.TemporaryDirectory(prefix=".chia-work-", dir=target_path)
        work = Path(temporary.name)
    try:
        promoters = work / "TSSgenesbed"
        chiageneTSS(genes, promoters, panther)
        combined = work / "linkschia"
        with combined.open("w", encoding="utf-8") as all_links:
            for sample, label, mango in samples:
                sample_dir = work / sample
                sample_dir.mkdir(parents=True, exist_ok=True)
                anchors = sample_dir / "anchors.bed"
                enhancer_overlap = sample_dir / "enhancers.intersect"
                gene_overlap = sample_dir / "genes.intersect"
                links = sample_dir / "links"
                celltypev2(mango, label, anchors)
                bedtoolIntersect(repo_root / "CREbedDBenhancers_10092018", anchors, enhancer_overlap)
                bedtoolIntersect(promoters, anchors, gene_overlap)
                chiagenesenhancers(enhancer_overlap, gene_overlap, links, panther, tissues)
                with links.open(encoding="utf-8") as source:
                    shutil.copyfileobj(source, all_links)
                if not keep_intermediates:
                    shutil.rmtree(sample_dir)
                print(f"ChIA-PET sample {sample} done", flush=True)
        final = work / "linksDBchia.complete"
        orderlinks(combined, final)
        os.replace(final, target_path / "linksDBchia")
        with (target_path / "linksDBchia").open(encoding="utf-8") as result:
            link_count = sum(1 for _ in result)
        print(f"ChIA-PET linksDBchia done: {link_count} links", flush=True)
        if link_count == 0:
            print("WARNING: ChIA-PET completed with zero enhancer-gene links; verify input data and mapping coverage.", file=sys.stderr, flush=True)
    finally:
        if temporary is not None:
            temporary.cleanup()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_path", type=Path, help="ChIA output directory, e.g. results/chia_pet")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--manifest", type=Path, help="sample TSV; default: chia_pet/samples.tsv")
    parser.add_argument("--genes", type=Path, help="Ensembl six-column table; default: chia_pet/hg19genes.txt")
    parser.add_argument("--check-inputs", action="store_true", help="validate inputs without running BEDTools")
    parser.add_argument("--keep-intermediates", action="store_true", help="retain anchors, overlaps and per-sample links")
    args = parser.parse_args(argv)
    repo = args.repo_root.resolve()
    manifest = args.manifest.resolve() if args.manifest else repo / "chia_pet" / "samples.tsv"
    genes = args.genes.resolve() if args.genes else repo / "chia_pet" / "hg19genes.txt"
    try:
        if args.check_inputs:
            _, _, samples = validate_inputs(repo, manifest, genes)
            print(f"ChIA-PET inputs OK: {len(samples)} sample(s)")
        else:
            run_workflow(args.target_path.resolve(), repo, manifest, genes, args.keep_intermediates)
    except (OSError, ValueError, subprocess.CalledProcessError, sqlite3.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
