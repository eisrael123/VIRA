"""Step 3: collapse clustered CAGE peaks down to main peak(s) per gene."""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class _Peak:
    line: str
    peak_pos: int
    depth: int


def choose_main_peaks(
    matched_bed: str,
    output_dir: str,
    cluster_distance: int = 100,
    distinct_distance: int = 200,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"Output_main_PEAKS_{os.path.basename(matched_bed)}.txt")

    genes: dict[str, list[_Peak]] = {}
    with open(matched_bed) as fh:
        for line in fh:
            line = line.rstrip("\n")
            cols = line.split("\t")
            peak_info = cols[3].split(";")
            peak_pos = int(peak_info[1])
            depth = int(cols[4])
            gene = cols[6]
            genes.setdefault(gene, []).append(_Peak(line=line, peak_pos=peak_pos, depth=depth))

    final_peaks: list[str] = []
    for gene in sorted(genes):
        peaks = sorted(genes[gene], key=lambda p: p.peak_pos)
        clusters: list[list[_Peak]] = []
        current = [peaks[0]]

        for peak in peaks[1:]:
            if abs(peak.peak_pos - current[-1].peak_pos) <= cluster_distance:
                current.append(peak)
            else:
                clusters.append(current)
                current = [peak]
        clusters.append(current)

        kept: list[_Peak] = []
        last_rep_pos: int | None = None
        for cluster in clusters:
            max_depth = max(p.depth for p in cluster)
            strongest = [p for p in cluster if p.depth == max_depth]
            rep_pos = strongest[0].peak_pos
            if last_rep_pos is None or abs(rep_pos - last_rep_pos) >= distinct_distance:
                kept.extend(strongest)
                last_rep_pos = rep_pos

        final_peaks.extend(p.line for p in kept)

    def sort_key(line: str) -> tuple[str, int]:
        cols = line.split("\t")
        return (cols[0], int(cols[1]))

    final_peaks.sort(key=sort_key)

    with open(output_file, "w") as out:
        for line in final_peaks:
            out.write(line + "\n")

    return output_file
