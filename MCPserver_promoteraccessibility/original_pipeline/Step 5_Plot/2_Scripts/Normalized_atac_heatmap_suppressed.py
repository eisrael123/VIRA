
import sys
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
import random
from scipy.ndimage import gaussian_filter1d

# ------------------
# INPUTS
# ------------------
bed_path = 'Suppressed_genes_DG75_BRRF1.txt_Output_main_PEAKS_Matched_TSS_100_TS_50_BRRF1_PEAKS_max_width_16.fraction_max_peaks_0.2.min_single_pos_depth_10.strand_MERGED.txt'
brrf1 = ['/Users/trangnguyen/Documents/8_Gene_expression_paper/Analysis/ATAC/4_Output_for_heatmap/1_insertion_bw_atac/BRRF_1-6_clean.ins.bw']
ctl = ['/Users/trangnguyen/Documents/8_Gene_expression_paper/Analysis/ATAC/4_Output_for_heatmap/1_insertion_bw_atac/Control_2-6_clean.ins.bw']

upstream_bases = 5000
downstream_bases = 5000
length = upstream_bases + downstream_bases

# ------------------
# LOAD REGIONS
# ------------------
regions = []
with open(bed_path) as bed_handle:
    for line in bed_handle:
        regions.append(line.strip().split('\t'))

random.shuffle(regions)
regions.sort(key=lambda x: float(x[5]))

# ------------------
# COVERAGE EXTRACTION
# ------------------
def extract_coverage(bigwig_paths, regions):
    m = np.zeros((len(regions), length))

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

            if strand == '-':
                vals = vals[::-1]

            m[ind] += vals

        bw.close()
        print(bw_path, "done")

    return m / len(bigwig_paths)

# ------------------
# EXTRACT DATA
# ------------------
ctl_coverage = extract_coverage(ctl, regions)
brrf1_coverage = extract_coverage(brrf1, regions)

# ------------------
# PI GLOBAL NORMALIZATION
# ------------------

total_ctl = np.sum(ctl_coverage)
total_brrf1 = np.sum(brrf1_coverage)

avg_total = (total_ctl + total_brrf1) / 2

ctl_coverage = (ctl_coverage / total_ctl) * avg_total
brrf1_coverage = (brrf1_coverage / total_brrf1) * avg_total

print("Normalization complete")
print("Total Control (before):", total_ctl)
print("Total BRRF1 (before):", total_brrf1)

# ------------------
# HEATMAPS
# ------------------
height = 8

for matrix, label in zip([ctl_coverage, brrf1_coverage], ["ctl", "BRRF1"]):
    plt.figure(figsize=(3, height))

    vmax_value = np.percentile(matrix, 99)

    plt.imshow(
        matrix,
        cmap='Reds',
        aspect='auto',
        interpolation='gaussian',
        vmax=vmax_value
    )

    plt.xticks([length // 2], ['TSS'])
    plt.yticks([])
    plt.tight_layout()
    plt.savefig(f'{label}_atac_dg75_{length}_BRRF1_test_SUPPRESSED.png', dpi=500)
    plt.close()

# ------------------
# SUMMATION PLOTS (SMOOTHED)
# ------------------

sigma = 50

ctl_sum = np.sum(ctl_coverage, axis=0)
brrf1_sum = np.sum(brrf1_coverage, axis=0)

ctl_smooth = gaussian_filter1d(ctl_sum, sigma=sigma)
brrf1_smooth = gaussian_filter1d(brrf1_sum, sigma=sigma)

ymax = max(ctl_smooth.max(), brrf1_smooth.max())

# ------------------
# COMBINED META-PLOT
# ------------------

plt.figure(figsize=(4, 3))

plt.plot(
    ctl_smooth,
    linewidth=1,
    label='Control'
)

plt.plot(
    brrf1_smooth,
    linewidth=1,
    label='BRRF1'
)

plt.xlabel("Position (bp)")
plt.ylabel("Normalized summed signal")
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig('summation_curve_combined_atac_dg75_BRRF1_test_SUPPRESSED.png', dpi=500)
plt.close()