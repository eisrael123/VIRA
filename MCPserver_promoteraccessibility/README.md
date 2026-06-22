# Promoter-accessibility ATAC pipeline (MCP server)

A reusable, path-free port of the original Step 2-5 scripts. Feed the inputs in,
get every intermediate file plus the nine Step-5 figures out. Steps 2-4 (the
Perl scripts) have been ported to Python; the three near-identical Step-5
plotting scripts have been consolidated into one.

## What it does

| Step | Module | In | Out |
|------|--------|----|-----|
| 2 | `pipeline/step2_tss_peaks.py` | peak BED + annotation BED + cutoffs | `Matched_*` / `Unmatched_*` |
| 3 | `pipeline/step3_main_peaks.py` | matched peaks + distances | `Output_main_PEAKS_*` |
| 4 | `pipeline/step4_extract_tpm.py` | DESeq2 `.xlsx` + main peaks | `*_induced/_suppressed/_nonsig.txt` (6-col plot inputs) + gene lists |
| 5 | `pipeline/step5_plot.py` | plot inputs + insertion bigwigs | 9 PNG figures |

`pipeline/run_pipeline.py` orchestrates all of it; `server.py` exposes it as a
single MCP tool, `run_promoter_accessibility_pipeline`.

## Install

```bash
cd "/Users/mac15/ethan/VIRA/MCPserver_promoteraccessibility"
uv sync
```

## Required inputs

- `--peaks-bed`: Step-1 peak-calling BED (e.g. the BRRF1 `*.strand_MERGED.bed`)
- `--annotation-bed`: `Step 2_ Extracting_TSS_peaks/Annotation file/Hg38_refgene.bed`
- `--deseq2-xlsx`: `Step 4_Extracting_TPM/Deseq2/DG75_BRRF1_test_vs_cntl_deseq2_results_genes.xlsx`
- `--test-bigwig` / `--ctl-bigwig`: insertion `.ins.bw` files (only needed for Step 5; on NGS4)

## Run from the command line

Steps 2-4 only (no bigwigs needed):

```bash
uv run python -m pipeline.run_pipeline \
  --peaks-bed "Step 1_Peak-caling/Peak_calling_output_files/BRRF1_PEAKS_max_width_16.fraction_max_peaks_0.2.min_single_pos_depth_10.strand_MERGED.bed" \
  --annotation-bed "Step 2_ Extracting_TSS_peaks/Annotation file/Hg38_refgene.bed" \
  --deseq2-xlsx "Step 4_Extracting_TPM/Deseq2/DG75_BRRF1_test_vs_cntl_deseq2_results_genes.xlsx" \
  --output-dir "pipeline_output" \
  --sample-label DG75_BRRF1 \
  --skip-plots
```

Full pipeline including figures (add the bigwigs and drop `--skip-plots`):

```bash
uv run python -m pipeline.run_pipeline \
  --peaks-bed "Step 1_Peak-caling/Peak_calling_output_files/BRRF1_PEAKS_max_width_16.fraction_max_peaks_0.2.min_single_pos_depth_10.strand_MERGED.bed" \
  --annotation-bed "Step 2_ Extracting_TSS_peaks/Annotation file/Hg38_refgene.bed" \
  --deseq2-xlsx "Step 4_Extracting_TPM/Deseq2/DG75_BRRF1_test_vs_cntl_deseq2_results_genes.xlsx" \
  --output-dir "pipeline_output" \
  --sample-label DG75_BRRF1 \
  --test-bigwig /path/to/BRRF_1-6_clean.ins.bw \
  --ctl-bigwig /path/to/Control_2-6_clean.ins.bw
```

Outputs land under `pipeline_output/step2`, `step3`, `step4`, and
`step5/figures` (9 PNGs).

## Run as an MCP server

```bash
uv run mcp dev server.py     # open the MCP Inspector to call the tool
uv run python server.py      # stdio server for an MCP client
```

Tunable parameters (defaults match the original scripts): `tss_cutoff=100`,
`ts_cutoff=50`, `cluster_distance=100`, `distinct_distance=200`, `log2fc=0.5`,
`padj=0.005`, `min_cntl_tpm=3`, `upstream_bases=5000`, `downstream_bases=5000`,
`sigma=50`.

## Assumptions

- Only one (BRRF1) peak file flows through to the figures.
- DESeq2 columns use `_cntl<N>_` / `_test<N>_` naming, with `log2FoldChange`
  and `padj` headers.
- The original `Step 1`-`Step 5` folders are kept untouched as reference data.
