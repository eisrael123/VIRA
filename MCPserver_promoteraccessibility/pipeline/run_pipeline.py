"""Run the promoter-accessibility pipeline from Step 2 through Step 5."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field

from pipeline.step2_tss_peaks import extract_tss_peaks
from pipeline.step3_main_peaks import choose_main_peaks
from pipeline.step4_extract_tpm import extract_tpm


@dataclass
class PipelineResult:
    matched_peaks: str
    unmatched_peaks: str
    main_peaks: str
    induced_plot_input: str
    suppressed_plot_input: str
    nonsig_plot_input: str
    induced_genes: str
    suppressed_genes: str
    nonsig_genes: str
    deseq2_tsv: str
    figures: list[str] = field(default_factory=list)

    def all_outputs(self) -> list[str]:
        return [
            self.matched_peaks,
            self.unmatched_peaks,
            self.main_peaks,
            self.induced_plot_input,
            self.suppressed_plot_input,
            self.nonsig_plot_input,
            self.induced_genes,
            self.suppressed_genes,
            self.nonsig_genes,
            self.deseq2_tsv,
            *self.figures,
        ]


def run_pipeline(
    peaks_bed: str,
    annotation_bed: str,
    deseq2_xlsx: str,
    output_dir: str,
    test_bigwigs: list[str] | None = None,
    ctl_bigwigs: list[str] | None = None,
    *,
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
) -> PipelineResult:
    for label, path in [
        ("peaks_bed", peaks_bed),
        ("annotation_bed", annotation_bed),
        ("deseq2_xlsx", deseq2_xlsx),
    ]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{label} not found: {path}")

    step2_dir = os.path.join(output_dir, "step2")
    step3_dir = os.path.join(output_dir, "step3")
    step4_dir = os.path.join(output_dir, "step4")
    figures_dir = os.path.join(output_dir, "step5", "figures")

    print("[Step 2] Extracting TSS peaks...")
    step2 = extract_tss_peaks(
        peak_bed=peaks_bed,
        annotation_bed=annotation_bed,
        output_dir=step2_dir,
        tss_cutoff=tss_cutoff,
        ts_cutoff=ts_cutoff,
    )

    print("[Step 3] Choosing main peaks...")
    main_peaks = choose_main_peaks(
        matched_bed=step2.matched_path,
        output_dir=step3_dir,
        cluster_distance=cluster_distance,
        distinct_distance=distinct_distance,
    )

    print("[Step 4] Extracting TPM and classifying genes...")
    step4 = extract_tpm(
        deseq2_xlsx=deseq2_xlsx,
        main_peaks=main_peaks,
        output_dir=step4_dir,
        log2fc=log2fc,
        padj=padj,
        min_cntl_tpm=min_cntl_tpm,
        sample_label=sample_label,
    )

    figures: list[str] = []
    if skip_plots:
        print("[Step 5] Skipped (skip_plots=True).")
    elif not test_bigwigs or not ctl_bigwigs:
        print("[Step 5] Skipped: no test/control bigwig files provided.")
    else:
        from pipeline.step5_plot import plot_class

        resolved_sample = sample_label or "sample"
        classes = [
            ("INDUCED", step4.induced_plot),
            ("SUPPRESSED", step4.suppressed_plot),
            ("NONSIG", step4.nonsig_plot),
        ]
        for class_label, bed_path in classes:
            print(f"[Step 5] Plotting {class_label}...")
            figures.extend(
                plot_class(
                    bed_path=bed_path,
                    test_bigwigs=test_bigwigs,
                    ctl_bigwigs=ctl_bigwigs,
                    output_dir=figures_dir,
                    class_label=class_label,
                    sample_label=resolved_sample,
                    test_label=test_label,
                    ctl_label=ctl_label,
                    upstream_bases=upstream_bases,
                    downstream_bases=downstream_bases,
                    sigma=sigma,
                )
            )

    return PipelineResult(
        matched_peaks=step2.matched_path,
        unmatched_peaks=step2.unmatched_path,
        main_peaks=main_peaks,
        induced_plot_input=step4.induced_plot,
        suppressed_plot_input=step4.suppressed_plot,
        nonsig_plot_input=step4.nonsig_plot,
        induced_genes=step4.induced_genes,
        suppressed_genes=step4.suppressed_genes,
        nonsig_genes=step4.nonsig_genes,
        deseq2_tsv=step4.tsv_path,
        figures=figures,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the promoter-accessibility pipeline.")
    parser.add_argument("--peaks-bed", required=True, help="Step-1 CAGE peak-calling BED file.")
    parser.add_argument("--annotation-bed", required=True, help="Gene annotation BED.")
    parser.add_argument("--deseq2-xlsx", required=True, help="DESeq2 results .xlsx file.")
    parser.add_argument("--output-dir", required=True, help="Directory for all outputs.")
    parser.add_argument("--test-bigwig", action="append", default=[], help="Test/BRRF1 insertion bigwig.")
    parser.add_argument("--ctl-bigwig", action="append", default=[], help="Control insertion bigwig.")
    parser.add_argument("--tss-cutoff", type=int, default=100)
    parser.add_argument("--ts-cutoff", type=int, default=50)
    parser.add_argument("--cluster-distance", type=int, default=100)
    parser.add_argument("--distinct-distance", type=int, default=200)
    parser.add_argument("--log2fc", type=float, default=0.5)
    parser.add_argument("--padj", type=float, default=0.005)
    parser.add_argument("--min-cntl-tpm", type=float, default=3.0)
    parser.add_argument("--upstream-bases", type=int, default=5000)
    parser.add_argument("--downstream-bases", type=int, default=5000)
    parser.add_argument("--sigma", type=int, default=50)
    parser.add_argument("--sample-label", default=None)
    parser.add_argument("--test-label", default="BRRF1")
    parser.add_argument("--ctl-label", default="ctl")
    parser.add_argument("--skip-plots", action="store_true")
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    result = run_pipeline(
        peaks_bed=args.peaks_bed,
        annotation_bed=args.annotation_bed,
        deseq2_xlsx=args.deseq2_xlsx,
        output_dir=args.output_dir,
        test_bigwigs=args.test_bigwig,
        ctl_bigwigs=args.ctl_bigwig,
        tss_cutoff=args.tss_cutoff,
        ts_cutoff=args.ts_cutoff,
        cluster_distance=args.cluster_distance,
        distinct_distance=args.distinct_distance,
        log2fc=args.log2fc,
        padj=args.padj,
        min_cntl_tpm=args.min_cntl_tpm,
        upstream_bases=args.upstream_bases,
        downstream_bases=args.downstream_bases,
        sigma=args.sigma,
        sample_label=args.sample_label,
        test_label=args.test_label,
        ctl_label=args.ctl_label,
        skip_plots=args.skip_plots,
    )

    print("\nDone. Generated files:")
    for path in result.all_outputs():
        print(f"  {path}")


if __name__ == "__main__":
    main()
