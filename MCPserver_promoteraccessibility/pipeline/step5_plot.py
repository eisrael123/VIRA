"""Step 5: ATAC insertion heatmaps and summation meta-plots."""

from __future__ import annotations

import os
import random

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyBigWig  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402


def _load_regions(bed_path: str) -> list[list[str]]:
    regions: list[list[str]] = []
    with open(bed_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped:
                regions.append(stripped.split("\t"))
    random.shuffle(regions)
    regions.sort(key=lambda x: float(x[5]))
    return regions


def _extract_coverage(
    bigwig_paths: list[str],
    regions: list[list[str]],
    upstream_bases: int,
    downstream_bases: int,
) -> np.ndarray:
    length = upstream_bases + downstream_bases
    matrix = np.zeros((len(regions), length))

    for bw_path in bigwig_paths:
        bw = pyBigWig.open(bw_path)
        for ind, region in enumerate(regions):
            chrom = region[0]
            start = int(float(region[1]))
            stop = int(float(region[2]))
            strand = region[3]

            middle = (start + stop) // 2
            start = middle - upstream_bases
            stop = middle + downstream_bases

            vals = bw.values(chrom, start, stop)
            vals = np.nan_to_num(np.abs(vals), nan=0)
            if strand == "-":
                vals = vals[::-1]
            matrix[ind] += vals
        bw.close()
        print(bw_path, "done")

    return matrix / len(bigwig_paths)


def plot_class(
    bed_path: str,
    test_bigwigs: list[str],
    ctl_bigwigs: list[str],
    output_dir: str,
    class_label: str,
    sample_label: str = "sample",
    test_label: str = "BRRF1",
    ctl_label: str = "ctl",
    upstream_bases: int = 5000,
    downstream_bases: int = 5000,
    sigma: int = 50,
    heatmap_height: int = 8,
) -> list[str]:
    os.makedirs(output_dir, exist_ok=True)
    length = upstream_bases + downstream_bases
    regions = _load_regions(bed_path)
    if not regions:
        print(f"No regions in {bed_path}; skipping {class_label}.")
        return []

    ctl_coverage = _extract_coverage(ctl_bigwigs, regions, upstream_bases, downstream_bases)
    test_coverage = _extract_coverage(test_bigwigs, regions, upstream_bases, downstream_bases)

    total_ctl = np.sum(ctl_coverage)
    total_test = np.sum(test_coverage)
    avg_total = (total_ctl + total_test) / 2
    ctl_coverage = (ctl_coverage / total_ctl) * avg_total
    test_coverage = (test_coverage / total_test) * avg_total

    print("Normalization complete")
    print("Total Control (before):", total_ctl)
    print(f"Total {test_label} (before):", total_test)

    written: list[str] = []
    for matrix, label in zip([ctl_coverage, test_coverage], [ctl_label, test_label]):
        plt.figure(figsize=(3, heatmap_height))
        vmax_value = np.percentile(matrix, 99)
        plt.imshow(matrix, cmap="Reds", aspect="auto", interpolation="gaussian", vmax=vmax_value)
        plt.xticks([length // 2], ["TSS"])
        plt.yticks([])
        plt.tight_layout()
        out_path = os.path.join(output_dir, f"{label}_atac_{sample_label}_{length}_{class_label}.png")
        plt.savefig(out_path, dpi=500)
        plt.close()
        written.append(out_path)

    ctl_smooth = gaussian_filter1d(np.sum(ctl_coverage, axis=0), sigma=sigma)
    test_smooth = gaussian_filter1d(np.sum(test_coverage, axis=0), sigma=sigma)

    plt.figure(figsize=(4, 3))
    plt.plot(ctl_smooth, linewidth=1, label="Control")
    plt.plot(test_smooth, linewidth=1, label=test_label)
    plt.xlabel("Position (bp)")
    plt.ylabel("Normalized summed signal")
    plt.legend(frameon=False)
    plt.tight_layout()
    summation_path = os.path.join(
        output_dir, f"summation_curve_combined_atac_{sample_label}_{class_label}.png"
    )
    plt.savefig(summation_path, dpi=500)
    plt.close()
    written.append(summation_path)

    return written
