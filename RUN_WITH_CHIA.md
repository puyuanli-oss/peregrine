# Run PEREGRINE with ChIA-PET

This small overlay adds a launcher, preflight checker and reproducible tests to
your **full** `chia-pet-compatibility-review` checkout. It does not replace your
scientific algorithms, annotations, dependencies or input data. The uploaded
`peregrine_review_small.zip` contains selected files; it is not a standalone
distribution of the full pipeline.

## 1. Install the overlay

Download `PEREGRINE_Run_Check.zip`. In your existing `peregrine` repository root:

```bash
unzip "$HOME/Downloads/PEREGRINE_Run_Check.zip" -d "$PWD"
```

The archive contains only new helper files, this guide and verification files.
It does not overwrite `steps.sh`, the ChIA code, or any annotation/input file.

## 2. Test the supplied ChIA implementation

Use Python 3.9 or newer and a real `bedtools` executable on PATH:

```bash
python3 verify_peregrine.py --scope chia --self-test
```

This creates and deletes a private synthetic dataset. It runs your actual
ChIA script with real BEDTools, checks eight expected ChIA links, checks that
ChIA evidence changes TAD selection from two rows to four, and checks final
TAD sorting when `tad/orderlinks.py` is present. No real MANGO data or production
manifest is required for this synthetic test.

For the full regression suite, including comparisons to the supplied original
Perl scripts:

```bash
python3 -m unittest discover -s verification_tests -v
```

Perl is needed for the original-code comparisons; BEDTools is needed for actual
intersection tests. Missing tools cause those tests to be skipped. A skipped
intersection test does not demonstrate that BEDTools ran.

## 3. Supply production inputs

Keep the shared references and standard eQTL/TAD inputs from your full checkout.
Use your existing verified `venv/bin/python` for the full pipeline. The full
preflight imports NumPy, pandas, pybedtools and pysam; it does not install them.

ChIA additionally requires:

- `chia_pet/hg19genes.txt`, which you have already generated and committed.
- Real processed, significant MANGO interaction files. Each must be plain text,
  with eight tab-separated columns: chromosome1, start1, end1, chromosome2,
  start2, end2, PET count and p-value. URLs, FASTQ files, gzip files and generic
  downloadable interaction BED files are not accepted as substitutes.
- `chia_pet/samples.tsv`, identifying those real files.

The manifest has exactly this tab-separated header:

```text
sample_id	cell_type	target_protein	mango_file
```

Add one row for each actual processed experiment. Relative paths in `mango_file`
are relative to the manifest directory. For example, a **real** file stored at
`chia_pet/data/HCT116_POLR2A.interactions.fdr.mango` can be listed as:

```text
HCT116_POLR2A	HCT116	POLR2A	data/HCT116_POLR2A.interactions.fdr.mango
```

Use the actual cell type and target for every input. Do not point several rows
at one file merely to stand in for missing experiments. The kit deliberately
does not install a pretend production manifest or production data.

To avoid regenerating raw sequencing data, obtain the lab's original processed
MANGO outputs. Otherwise, the external MANGO workflow must generate them before
PEREGRINE can consume them. A manifest describes files; it does not download or
generate them. FASTQ and MANGO processing were upstream of Caitlin's Perl code.

## 4. Check readiness, then run

From the full repository root, with the working environment and BEDTools loaded:

```bash
bash run_peregrine.sh results_with_chia --check-inputs
```

This checks source syntax, imports, BEDTools, expected standard input paths and
the actual ChIA manifest/gene/MANGO contents. It exits without generating links.
The standard input checks confirm presence, not every possible content error in
the older eQTL/hierarchical/TAD implementations.

After the check passes, create a fresh output directory and run:

```bash
mkdir -p runs
PEREGRINE_OUT=$(mktemp -d "$PWD/runs/peregrine_XXXXXXXX")
printf 'Output directory: %s\n' "$PEREGRINE_OUT"
bash run_peregrine.sh "$PEREGRINE_OUT"
```

Use a new directory for a full run: older TAD stages may append intermediate
files on a rerun. For a compute cluster, put the final launcher command inside
your batch script with the same environment and resources you already use.

The new launcher always includes ChIA. Its defaults are
`chia_pet/samples.tsv` and `chia_pet/hg19genes.txt`; override them with
`--chia-manifest FILE` and `--chia-genes FILE` when needed.

The existing `steps.sh` is unchanged. Its direct equivalent remains:

```bash
bash steps.sh results_with_chia --chia-manifest chia_pet/samples.tsv
```

Calling `steps.sh` without `--chia-manifest` skips ChIA. It also requires an
existing `venv`; the old root README's claim that it creates one is outdated.

Execution order is eQTL, hierarchical TAD, ChIA-PET, then linear TAD. In linear
TAD, `cutdowntad.py` uses the union of ChIA, eQTL and hierarchical evidence;
`orderlinks.py` writes the sorted final TAD links. ChIA output is
`OUTPUT/chia_pet/linksDBchia` with four headerless fields:
enhancer ID, PANTHER gene ID, tissue code and assay `1`.

## Update an existing standard run

If standard stages have already completed in an output directory, supply that
same directory to both commands:

```bash
bash run_peregrine.sh OUTPUT_DIRECTORY --only-chia
bash run_peregrine.sh OUTPUT_DIRECTORY --final-only
```

The first produces ChIA links; the second repeats final TAD filtering and
sorting with those links. It requires the existing eQTL, hierarchical and TAD
candidate outputs. It does not rerun the earlier standard stages. Both commands
stop on failure. A successful run with no supported ChIA links is possible and
is reported explicitly by the ChIA code.

## What was verified here

`verification_tests/validation.log` and `validation_summary.json` record the
checks against the uploaded source. Fixtures are synthetic. Original Perl is
executed for byte comparisons; wrapper dispatch tests use recorded calls and
are not biological validation. Actual BEDTools tests were unavailable in the
review environment and are provided for execution on your machine.

Passing controlled tests supports compatibility of the tested linking rules.
It does not establish identical full historical results: those also depend on
the exact MANGO inputs, parameters and original gene annotation/order.
