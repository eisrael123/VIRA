"""MCP server exposing the promoter-accessibility pipeline as background jobs.

Run with:
    uv run mcp dev server.py      # interactive inspector
    uv run server.py              # stdio server (for an MCP client)
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
import uuid
from pathlib import Path

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("promoter-accessibility")

_JOBS: dict[str, "PipelineJob"] = {}


class PipelineJob:
    def __init__(
        self,
        job_id: str,
        pid: int,
        process: subprocess.Popen,
        command: list[str],
        output_dir: str,
        log_path: str,
        started_at: float,
    ) -> None:
        self.job_id = job_id
        self.pid = pid
        self.process = process
        self.command = command
        self.output_dir = output_dir
        self.log_path = log_path
        self.started_at = started_at


@mcp.tool()
def start_step2_tss_peaks(
    peak_bed: str,
    annotation_bed: str,
    output_dir: str,
    tss_cutoff: int = 100,
    ts_cutoff: int = 50,
) -> dict:
    """Start Step 2 in the background: match CAGE peaks to annotated TSS windows."""
    missing = _missing_paths([peak_bed, annotation_bed])
    if missing:
        return {"started": False, "error": "Required input file(s) are missing.", "missing_files": missing}

    command = [
        sys.executable,
        "-m",
        "pipeline.run_step2",
        "--peak-bed",
        peak_bed,
        "--annotation-bed",
        annotation_bed,
        "--output-dir",
        output_dir,
        "--tss-cutoff",
        str(tss_cutoff),
        "--ts-cutoff",
        str(ts_cutoff),
    ]
    result = _launch_background_job(command=command, output_dir=output_dir, job_kind="step2_tss_peaks")
    result["next_step"] = "After completion, read the log for matched_path and pass it to start_step3_main_peaks(...)."
    return result


@mcp.tool()
def start_step3_main_peaks(
    matched_bed: str,
    output_dir: str,
    cluster_distance: int = 100,
    distinct_distance: int = 200,
) -> dict:
    """Start Step 3 in the background: choose main CAGE peak(s) from Step-2 output."""
    missing = _missing_paths([matched_bed])
    if missing:
        return {"started": False, "error": "Required input file(s) are missing.", "missing_files": missing}

    command = [
        sys.executable,
        "-m",
        "pipeline.run_step3",
        "--matched-bed",
        matched_bed,
        "--output-dir",
        output_dir,
        "--cluster-distance",
        str(cluster_distance),
        "--distinct-distance",
        str(distinct_distance),
    ]
    result = _launch_background_job(command=command, output_dir=output_dir, job_kind="step3_main_peaks")
    result["next_step"] = "After completion, read the log for main_peaks and pass it to start_step4_extract_tpm(...)."
    return result


@mcp.tool()
def start_step4_extract_tpm(
    deseq2_xlsx: str,
    main_peaks: str,
    output_dir: str,
    log2fc: float = 0.5,
    padj: float = 0.005,
    min_cntl_tpm: float = 3.0,
    sample_label: str | None = None,
) -> dict:
    """Start Step 4 in the background: classify DESeq2 genes and create Step-5 inputs."""
    missing = _missing_paths([deseq2_xlsx, main_peaks])
    if missing:
        return {"started": False, "error": "Required input file(s) are missing.", "missing_files": missing}

    command = [
        sys.executable,
        "-m",
        "pipeline.run_step4",
        "--deseq2-xlsx",
        deseq2_xlsx,
        "--main-peaks",
        main_peaks,
        "--output-dir",
        output_dir,
        "--log2fc",
        str(log2fc),
        "--padj",
        str(padj),
        "--min-cntl-tpm",
        str(min_cntl_tpm),
    ]
    if sample_label:
        command.extend(["--sample-label", sample_label])
    result = _launch_background_job(command=command, output_dir=output_dir, job_kind="step4_extract_tpm")
    result["next_step"] = "After completion, read the log for the three plot input files and pass them to start_step5_plot(...)."
    return result


@mcp.tool()
def start_step5_plot(
    induced_bed: str,
    suppressed_bed: str,
    nonsig_bed: str,
    output_dir: str,
    test_bigwigs: list[str],
    ctl_bigwigs: list[str],
    sample_label: str = "sample",
    test_label: str = "BRRF1",
    ctl_label: str = "ctl",
    upstream_bases: int = 5000,
    downstream_bases: int = 5000,
    sigma: int = 50,
) -> dict:
    """Start only Step 5 plotting in the background from Step-4 plot inputs."""
    missing = _missing_paths([induced_bed, suppressed_bed, nonsig_bed, *(test_bigwigs or []), *(ctl_bigwigs or [])])
    if not test_bigwigs:
        missing.append("<missing test_bigwigs>")
    if not ctl_bigwigs:
        missing.append("<missing ctl_bigwigs>")
    if missing:
        return {"started": False, "error": "Required input file(s) are missing.", "missing_files": missing}

    command = [
        sys.executable,
        "-m",
        "pipeline.run_step5",
        "--induced-bed",
        induced_bed,
        "--suppressed-bed",
        suppressed_bed,
        "--nonsig-bed",
        nonsig_bed,
        "--output-dir",
        output_dir,
        "--sample-label",
        sample_label,
        "--test-label",
        test_label,
        "--ctl-label",
        ctl_label,
        "--upstream-bases",
        str(upstream_bases),
        "--downstream-bases",
        str(downstream_bases),
        "--sigma",
        str(sigma),
    ]
    for bw in test_bigwigs:
        command.extend(["--test-bigwig", bw])
    for bw in ctl_bigwigs:
        command.extend(["--ctl-bigwig", bw])

    return _launch_background_job(command=command, output_dir=output_dir, job_kind="step5_plot")


@mcp.tool()
def check_pipeline_status(job_id: str) -> dict:
    """Check whether a background pipeline job is still running."""
    job = _JOBS.get(job_id)
    if job is None:
        return {
            "found": False,
            "error": "Unknown job id. This can happen if the MCP server restarted.",
            "job_id": job_id,
        }

    return_code = job.process.poll()
    if return_code is None:
        status = "running"
    elif return_code == 0:
        status = "completed"
    else:
        status = "failed"

    return {
        "found": True,
        "job_id": job_id,
        "pid": job.pid,
        "status": status,
        "return_code": return_code,
        "elapsed_seconds": round(time.time() - job.started_at, 1),
        "output_dir": job.output_dir,
        "log_path": job.log_path,
    }


@mcp.tool()
def read_pipeline_log(job_id: str, tail_lines: int = 80) -> dict:
    """Read the tail of a background job log."""
    job = _JOBS.get(job_id)
    if job is None:
        return {
            "found": False,
            "error": "Unknown job id. This can happen if the MCP server restarted.",
            "job_id": job_id,
        }
    return {
        "found": True,
        "job_id": job_id,
        "log_path": job.log_path,
        "tail": _tail(job.log_path, tail_lines),
    }


@mcp.tool()
def list_pipeline_outputs(output_dir: str) -> dict:
    """List generated pipeline outputs under an output directory."""
    root = Path(output_dir)
    if not root.exists():
        return {"found": False, "error": f"Output directory does not exist: {output_dir}"}

    files = sorted(str(path) for path in root.rglob("*") if path.is_file())
    figures = [path for path in files if path.lower().endswith(".png")]
    summation_curves = [
        path for path in figures if Path(path).name.startswith("summation_curve_combined_atac_")
    ]

    return {
        "found": True,
        "output_dir": output_dir,
        "file_count": len(files),
        "figure_count": len(figures),
        "summation_curve_count": len(summation_curves),
        "figures": figures,
        "summation_curves": summation_curves,
        "files": files,
    }


def _launch_background_job(command: list[str], output_dir: str, job_kind: str) -> dict:
    os.makedirs(output_dir, exist_ok=True)
    job_id = uuid.uuid4().hex[:12]
    job_dir = os.path.join(output_dir, ".pipeline_jobs")
    os.makedirs(job_dir, exist_ok=True)
    log_path = os.path.join(job_dir, f"{job_id}.{job_kind}.log")

    with open(log_path, "w") as log:
        log.write("$ " + _quote_command(command) + "\n\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=os.path.dirname(os.path.abspath(__file__)),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    _JOBS[job_id] = PipelineJob(
        job_id=job_id,
        pid=process.pid,
        process=process,
        command=command,
        output_dir=output_dir,
        log_path=log_path,
        started_at=time.time(),
    )

    return {
        "started": True,
        "job_id": job_id,
        "pid": process.pid,
        "status": "running",
        "output_dir": output_dir,
        "log_path": log_path,
        "command": _quote_command(command),
        "next_steps": [
            f"Call check_pipeline_status(job_id='{job_id}') until status is completed or failed.",
            f"Call read_pipeline_log(job_id='{job_id}') to inspect progress.",
            f"Call list_pipeline_outputs(output_dir='{output_dir}') after completion.",
        ],
    }


def _missing_paths(paths: list[str]) -> list[str]:
    return [path for path in paths if path.startswith("<missing") or not os.path.isfile(path)]


def _tail(path: str, lines: int) -> str:
    with open(path) as fh:
        all_lines = fh.readlines()
    return "".join(all_lines[-lines:])


def _quote_command(command: list[str]) -> str:
    import shlex

    return " ".join(shlex.quote(part) for part in command)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
