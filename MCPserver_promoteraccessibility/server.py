"""MCP server exposing the promoter-accessibility ATAC pipeline as one tool.

Run with:
    uv run mcp dev server.py      # interactive inspector
    uv run server.py              # stdio server (for an MCP client)
"""

from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from pipeline.run_pipeline import run_pipeline

mcp = FastMCP("promoter-accessibility")


@mcp.tool()
def run_promoter_accessibility_pipeline(
    peaks_bed: str,
    annotation_bed: str,
    deseq2_xlsx: str,
    output_dir: str,
    test_bigwigs: list[str] | None = None,
    ctl_bigwigs: list[str] | None = None,
    tss_cutoff: int = 100,
    ts_cutoff: int = 50,
    cluster_distance: int = 100,
    distinct_distance: int = 200,
    log2fc: float = 0.5,
    padj: float = 0.005,
    min_cntl_tpm: float = 3.0,
    upstream_bases: int = 5000,
    downstream_bases: int = 5000,
    sigma: int = 50,
    sample_label: str | None = None,
    test_label: str = "BRRF1",
    ctl_label: str = "ctl",
    skip_plots: bool = False,
) -> dict:
    """Run Steps 2-5 of the ATAC promoter-accessibility pipeline.

    Feed in the Step-1 peak BED, the gene annotation BED, the DESeq2 .xlsx, and
    (optionally) the test/control insertion bigwigs. The tool writes every
    intermediate file plus the nine Step-5 figures under ``output_dir`` and
    returns the paths to all generated files.

    Set ``skip_plots=True`` (or omit the bigwigs) to run Steps 2-4 only.
    """
    result = run_pipeline(
        peaks_bed=peaks_bed,
        annotation_bed=annotation_bed,
        deseq2_xlsx=deseq2_xlsx,
        output_dir=output_dir,
        test_bigwigs=test_bigwigs,
        ctl_bigwigs=ctl_bigwigs,
        tss_cutoff=tss_cutoff,
        ts_cutoff=ts_cutoff,
        cluster_distance=cluster_distance,
        distinct_distance=distinct_distance,
        log2fc=log2fc,
        padj=padj,
        min_cntl_tpm=min_cntl_tpm,
        upstream_bases=upstream_bases,
        downstream_bases=downstream_bases,
        sigma=sigma,
        sample_label=sample_label,
        test_label=test_label,
        ctl_label=ctl_label,
        skip_plots=skip_plots,
    )

    return {
        "matched_peaks": result.matched_peaks,
        "unmatched_peaks": result.unmatched_peaks,
        "main_peaks": result.main_peaks,
        "induced_plot_input": result.induced_plot_input,
        "suppressed_plot_input": result.suppressed_plot_input,
        "nonsig_plot_input": result.nonsig_plot_input,
        "induced_genes": result.induced_genes,
        "suppressed_genes": result.suppressed_genes,
        "nonsig_genes": result.nonsig_genes,
        "deseq2_tsv": result.deseq2_tsv,
        "figures": result.figures,
        "all_outputs": result.all_outputs(),
    }


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
