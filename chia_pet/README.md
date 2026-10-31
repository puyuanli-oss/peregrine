# ChIA-PET links compatible with PEREGRINE

This stage ports the original PEREGRINE Perl linking workflow to Python.
It starts with processed Mango interaction files. FASTQ alignment and Mango
processing remain upstream prerequisites, as in the original workflow.
No changes to requirements.txt or to the eQTL/hierarchical-TAD algorithms
are required. This stage uses Python 3.9+ and the existing BEDTools command.

## Inputs

Use hg19/GRCh37 consistently. This stage does not lift over coordinates.

Shared inputs at the repository root:

- CREbedDBenhancers_10092018: headerless BED4, chromosome/start/end/enhancer ID.
- pantherGeneList.txt: PANTHER long ID and versionless Ensembl gene ID in
  the first two tab-separated fields; additional columns are allowed.
- tissuetable_10092018.txt: tissue code and name in the first two tab-separated fields.

ChIA-specific inputs:

1. An Ensembl table at chia_pet/hg19genes.txt, or a path supplied with --genes.
   The first six tab-separated columns must be in this order:

   ~~~text
   Gene stable ID	Transcription start site (TSS)	Gene end (bp)	Gene start (bp)	Strand	Chromosome/scaffold name
   ENSG00000000001	1000	1500	1000	1	1
   ENSG00000000002	2000	2000	1800	-1	10
   ~~~

   The shown header is optional; trailing extra columns are ignored.
   Chromosome names in this table omit the chr prefix. Strands are 1/-1.
   Gene IDs are matched exactly; no version removal is added.
   The script keeps the last gene-table record for each mapped gene.
   A promoter with a negative start is rejected rather than clipped.

2. Headerless Mango interaction files with exactly eight tab-separated fields:

   ~~~text
   chr1	100	200	chr1	1000	1100	10	0.01
   ~~~

   Fields are chromosome1/start1/end1/chromosome2/start2/end2/PET count/p-value.
   Coordinates must be valid nonnegative BED intervals. PET counts are
   nonnegative integers; p-values are finite numbers from 0 through 1.
   No additional PET-count or significance cutoff is applied.
   Enhancers, generated promoters and Mango anchors must use matching chromosome names.

3. A sample manifest. Copy samples.example.tsv to samples.tsv and replace its
   example rows with actual experiments. It is a TSV with exactly this header:

   ~~~text
   sample_id	cell_type	target_protein	mango_file
   HCT116_POLR2A	HCT116	POLR2A	data/HCT116_POLR2A.interactions.fdr.mango
   ~~~

   Relative Mango paths are resolved against the manifest's directory, so this
   example reads chia_pet/data/HCT116_POLR2A.interactions.fdr.mango when the
   manifest is chia_pet/samples.tsv. Absolute paths and spaces in paths work.
   Sample IDs must be unique and use letters/digits/dots/underscores/hyphens,
   beginning with a letter or digit. Cell/protein metadata is explicit:
   it is not guessed from a filename.

   Cell matching is case-sensitive. Hyphens are removed from the cell name
   and from tissue-table names; for example HeLa-S3 matches HeLaS3.
   Cell names otherwise use letters/digits/dots. Protein labels cannot
   contain whitespace. Separate experiments are joined independently.

## Running

From the repository root, validate without generating links:

~~~bash
./venv/bin/python chia_pet/ChIA_PET.py results/chia_pet \
  --manifest chia_pet/samples.tsv \
  --genes chia_pet/hg19genes.txt \
  --check-inputs
~~~

Run only ChIA:

~~~bash
./venv/bin/python chia_pet/ChIA_PET.py results/chia_pet \
  --manifest chia_pet/samples.tsv \
  --genes chia_pet/hg19genes.txt
~~~

Run the complete pipeline including ChIA:

~~~bash
bash steps.sh results \
  --chia-manifest chia_pet/samples.tsv \
  --chia-genes chia_pet/hg19genes.txt
~~~

The order is eQTL, hierarchical TAD, ChIA-PET, then linear TAD.
ChIA inputs are checked before the expensive earlier stages run.

The existing command still runs the three existing stages:

~~~bash
bash steps.sh results
~~~

ChIA is enabled only when --chia-manifest is supplied. A previous linksDBchia
in the same output directory is never implicitly used by a disabled run.
Standalone ChIA defaults to chia_pet/samples.tsv and chia_pet/hg19genes.txt
under its repository root. Use --repo-root to select a different shared
input root. Explicit --manifest/--genes paths are relative to the caller's
working directory. Output paths may be absolute or contain spaces.

To retain the original-style intermediates, add --keep-intermediates to the
standalone ChIA command. They are written under results/chia_pet/intermediates.
By default, per-sample anchors and intersections are deleted after that
sample is linked; other temporary work is cleaned up at the end. Final sorting
and deduplication uses a temporary SQLite file on the output filesystem.
The final linksDBchia is replaced only after the whole ChIA run succeeds.
If a rerun fails, an older final file may still exist; it is not evidence of
success. The top-level launcher stops and does not pass it into TAD.

## Preserved linking logic

- Each Mango interaction becomes two BED5 anchors with a shared interaction name.
- Promoters are [TSS-600,TSS) on strand 1 and [TSS,TSS+600) on strand -1.
  The original arithmetic is preserved without an extra coordinate shift.
- Only genes in the PANTHER map get promoters.
- Both enhancer and promoter overlap use exactly:

  ~~~bash
  bedtools intersect -wa -wb -F 0.5 -f 0.5 -e
  ~~~

  At least 50% of EITHER feature or anchor must overlap, not necessarily both.
- A promoter on one anchor links to enhancers on its interacting PARTNER.
  Merely sharing the same anchor is insufficient unless the interaction is
  explicitly a self-interaction.
- Graphs are scoped by cell/protein metadata within each experiment.
- Final links are unique enhancer/PANTHER/tissue triples, sorted lexically
  by those fields. Assay is 1. There is no new p-value aggregation rule.

## Output and the original README discrepancy

The final file is results/chia_pet/linksDBchia:

~~~text
ENHANCER_ID	PANTHER_LONG_ID	TISSUE_CODE	1
~~~

It has FOUR tab-separated fields and no header. The original executable
chiagenesenhancers.pl and orderlinks.pl both write this format, although the
original README shows a five-field example containing a p-value. The rewrite
follows the executable code. Mango p-values remain in retained anchor/overlap
files when --keep-intermediates is used, but are not added to the final links.

## Linear TAD integration

The top-level launcher explicitly forwards this run's linksDBchia through:

~~~bash
bash tad/tad.sh results/tad --chia-links results/chia_pet/linksDBchia
~~~

Only cutdowntad.py consumes this extra argument:

~~~bash
./venv/bin/python tad/cutdowntad.py results/tad \
  --chia-links results/chia_pet/linksDBchia
~~~

A TAD row survives when its enhancer/PANTHER pair occurs in ChIA OR eQTL OR
hierarchical TAD. Evidence tissue need not match the TAD tissue. The TAD row
retains its own tissue and assay 3. An explicitly requested missing ChIA file
is an error; omitting the option retains the previous eQTL/hierarchy-only
behavior.

If TAD links were already generated without ChIA, run the new cutdowntad.py
command above, then tad/orderlinks.py with the same results/tad directory.
Adding ChIA support changes TAD filtering, so generating linksDBchia alone
does not update an old linksDBtad.
