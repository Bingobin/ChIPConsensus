# ChIP Consensus Peaks Toolkit

This program is a lightweight ChIP peak integration pipeline for multi-sample narrowPeak datasets. It takes peak calls from multiple samples (for example different patients, cell lines, or TF ChIP runs), merges overlapping regions across all samples into a consensus peak set, then quantifies sample-level support for each consensus region. The goal is to provide a reproducible, scriptable way to move from per-sample peak files to comparable cross-sample overlap results.

The pipeline is implemented in pure Python (no bedtools required), supports CSV/TSV config input with named sample columns (`ID`, `PeakFile`), and produces both tabular outputs and UpSet visualization inputs. It is suitable for downstream tasks such as overlap analysis, recurrent peak discovery, and sample-by-sample presence/absence pattern comparison.

Main outputs:

1. merged consensus BED (`chr`, `start`, `end`)
2. support matrix (`support_count` + one `0/1` column per sample)

## Input config format

Config must be CSV or TSV, with required columns `ID` and `PeakFile`. Columns are
matched by name, not position: their order can vary, and extra columns (such as
`Big`) are ignored. Column names are case-insensitive and surrounding whitespace
is ignored. The original two-column format remains supported.

CSV example:

```text
ID,PeakFile
SampleID1,/absolute/path/to/sample1.narrowPeak
SampleID2,/absolute/path/to/sample2.narrowPeak
```

Extra-column example (also valid):

```text
ID,PeakFile,Big
SampleID1,/absolute/path/to/sample1.narrowPeak,/absolute/path/to/sample1.bw
SampleID2,/absolute/path/to/sample2.narrowPeak,/absolute/path/to/sample2.bw
```

TSV example:

```text
ID	PeakFile
SampleID1	/absolute/path/to/sample1.narrowPeak
SampleID2	/absolute/path/to/sample2.narrowPeak
```

Lines starting with `#` and blank lines are ignored.
Each sample must have a non-empty, unique `ID` and a non-empty `PeakFile` pointing
to an existing file. Ignored columns do not need to contain existing file paths.

## Step-by-step scripts

Step scripts are in `src/`.

1. Validate config and summarize samples:

```bash
python3 src/step1_validate_config.py \
  --config samples.config.tsv \
  --out-manifest samples.manifest.tsv
```

2. Build merged consensus regions:

```bash
python3 src/step2_build_consensus.py \
  --config samples.config.tsv \
  --out-consensus consensus.bed
```

3. Build support matrix:

```bash
python3 src/step3_build_matrix.py \
  --config samples.config.tsv \
  --out-matrix consensus.support.tsv
```

4. Plot UpSet figure:

```bash
python3 src/step4_plot_upset.py \
  --matrix consensus.support.tsv \
  --out-pdf consensus.upset.pdf \
  --top-n 20
```

4R. Plot UpSet figure (R + ComplexHeatmap):

```bash
Rscript src/step4_plot_upset.R \
  --matrix consensus.support.tsv \
  --out-pdf consensus.upset.R.pdf \
  --top-n 20
```

## One-command pipeline

Main script is at project root: `chip_consensus.py`.

```bash
python3 chip_consensus.py \
  --config samples.config.tsv \
  --prefix TFs \
  --out-dir results
```

Outputs will be written to `results/`.
Generated files:
- `TFs.samples.manifest.tsv`
- `TFs.consensus.bed`
- `TFs.consensus.support.tsv`
- `TFs.consensus.upset.pdf` (default: auto-plot, Python version)
- `TFs.consensus.upset.R.pdf` (default: auto-plot, R version)

Notes for one-command pipeline:
- It prints step-by-step runtime logs (`Step 1/4` ... `Step 4/4`).
- By default it tries both Python and R UpSet plotting.
- If either UpSet plotting step fails (e.g. missing matplotlib or R package), pipeline continues and prints warning logs.
- To skip plotting: add `--no-upset`.

## Output matrix columns

`chr  start  end  support_count  SampleID1  SampleID2  ...`

Support rule: if a sample has any overlap (>=1 bp) with a consensus region, mark as `1`; otherwise `0`.
