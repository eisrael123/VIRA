"""
uv run mcp install server.py
"""

from mcp.server.fastmcp import FastMCP
import gget
import os
import re
from datetime import datetime

# Create an MCP server
mcp = FastMCP("gget-explorer", json_response=True)

DEFAULT_GGET_OUTPUT_DIR = "/Volumes/TUNGSACore3/VIRA_database/gget_data"

@mcp.resource("documentation://gget-virus")
def get_gget_virus_docs() -> str:
    """Provides documentation for running gget virus queries."""
    # You can point this to a local copy of the agent docs from the VirBench repo
    docs_path = os.path.join(os.path.dirname(__file__), "virus.md")
    
    if os.path.exists(docs_path):
        with open(docs_path, "r") as f:
            return f.read()
            
    return """
    ## gget virus Quick Reference
    Arguments:
    - virus: Target pathogen name, TaxID, accession, accession list, or accession file path
    - host: Host organism taxonomy name (e.g., 'human')
    - nuc_completeness: 'complete' or 'partial'
    - min_release_date: 'YYYY-MM-DD'
    - max_release_date: 'YYYY-MM-DD'
    - min_seq_length: integer min base pairs
    - max_seq_length: integer max base pairs
    
    Returns: Paths to FASTA, CSV, and JSONL files saved by gget virus.
    """

@mcp.tool()
def get_gget_virus_documentation() -> str:
    """Return the gget virus documentation."""
    return get_gget_virus_docs()


@mcp.tool()
def query_viral_sequences(
    virus: str,
    host: str = None,
    nuc_completeness: str = None,
    min_seq_length: int = None,
    max_seq_length: int = None,
    min_collection_date: str = None,
    max_collection_date: str = None,
    min_release_date: str = None,
    max_release_date: str = None,
    geographic_location: str = None,
    source_database: str = None,
    is_accession: bool = False,
    is_sars_cov2: bool = False,
    is_alphainfluenza: bool = False,
    outfolder: str = None,
) -> str:
    """
    Programmatically queries and filters the NCBI Virus database using gget virus.
    Always consult 'documentation://gget-virus' to understand structural parameters.

    Results are saved to an output folder as FASTA, CSV, and JSONL files.
    """
    try:
        if outfolder is None:
            safe_virus = re.sub(r"[^A-Za-z0-9_.-]+", "_", virus).strip("_") or "virus"
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            outfolder = os.path.join(
                DEFAULT_GGET_OUTPUT_DIR,
                f"{safe_virus}_{timestamp}",
            )

        result_path = os.path.abspath(outfolder)

        # Run gget virus via the Python API. gget writes output files to result_path.
        gget.virus(
            virus=virus,
            is_accession=is_accession,
            outfolder=result_path,
            host=host,
            min_seq_length=min_seq_length,
            max_seq_length=max_seq_length,
            nuc_completeness=nuc_completeness,
            geographic_location=geographic_location,
            min_collection_date=min_collection_date,
            max_collection_date=max_collection_date,
            min_release_date=min_release_date,
            max_release_date=max_release_date,
            source_database=source_database,
            is_sars_cov2=is_sars_cov2,
            is_alphainfluenza=is_alphainfluenza,
        )

        return (
            "gget virus query completed successfully.\n"
            f"Results were saved to: {result_path}\n"
            "Expected output files include FASTA, CSV metadata, and JSONL metadata when records are found."
        )

    except Exception as e:
        return f"Error executing gget virus query: {str(e)}"

@mcp.tool()
def read_gget_result_file(relative_path: str, max_chars: int = 10000) -> str:
    """Read a preview of a saved gget result file by path relative to the output directory."""
    base_path = os.path.abspath(DEFAULT_GGET_OUTPUT_DIR)
    file_path = os.path.abspath(os.path.join(base_path, relative_path))

    if not file_path.startswith(base_path):
        raise ValueError("File must be inside the gget output directory.")

    with open(file_path, "r") as f:
        return f.read(max_chars)

# Run with streamable HTTP transport
if __name__ == "__main__":
    mcp.run(transport="stdio")