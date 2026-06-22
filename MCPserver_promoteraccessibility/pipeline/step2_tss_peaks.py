"""Step 2: extract CAGE peaks that fall in a TSS window of annotated genes."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass


@dataclass
class Step2Result:
    matched_path: str
    unmatched_path: str


def _gene_base(gene_full: str) -> str:
    gene_full = re.sub(r"\s+", "", gene_full)
    return re.sub(r"[_-].*$", "", gene_full)


def extract_tss_peaks(
    peak_bed: str,
    annotation_bed: str,
    output_dir: str,
    tss_cutoff: int = 100,
    ts_cutoff: int = 50,
) -> Step2Result:
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.basename(peak_bed)

    matched_path = os.path.join(output_dir, f"Matched_TSS_{tss_cutoff}_TS_{ts_cutoff}_{filename}")
    unmatched_path = os.path.join(output_dir, f"Unmatched_TSS_{tss_cutoff}_TS_{ts_cutoff}_{filename}")

    annotations: list[list[str]] = []
    with open(annotation_bed) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 6:
                annotations.append(cols)

    with open(peak_bed) as fh, open(matched_path, "w") as out_matched, open(
        unmatched_path, "w"
    ) as out_unmatched:
        for line in fh:
            line = line.rstrip("\n")
            peak = line.split("\t")
            if len(peak) <= 4:
                continue

            chrom = peak[0]
            start_peak = int(peak[1])
            end_peak = int(peak[2])
            peak_strand = peak[5]

            matched = False
            seen_gene: set[str] = set()

            for anno in annotations:
                if chrom != anno[0] or peak_strand != anno[5]:
                    continue

                gene_base = _gene_base(anno[3])

                if anno[5] == "+":
                    upstream_tss = int(anno[1]) - tss_cutoff
                    thick_start = int(anno[6]) + ts_cutoff
                    in_window = start_peak > upstream_tss and end_peak < thick_start
                elif anno[5] == "-":
                    upstream_tss = int(anno[2]) + tss_cutoff
                    thick_start = int(anno[7]) - ts_cutoff
                    in_window = start_peak > thick_start and end_peak < upstream_tss
                else:
                    continue

                if in_window:
                    if gene_base in seen_gene:
                        continue
                    seen_gene.add(gene_base)
                    out_matched.write("\t".join(peak + [gene_base]) + "\n")
                    matched = True

            if not matched:
                out_unmatched.write("\t".join(peak) + "\n")

    return Step2Result(matched_path=matched_path, unmatched_path=unmatched_path)
