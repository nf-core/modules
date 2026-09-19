#!/usr/bin/env python3

import os
import platform

import pandas as pd

# Configure variables from Nextflow
prefix = "${task.ext.prefix}" if "${task.ext.prefix}" != "null" else "${meta.id}"
length_tool = "${length_tool}" if "${length_tool}" != "null" else "unknown"
length_tsv = "${length_tsv}"
content_tsv = "${content_tsv}" if "${content_tsv}" != "null" else ""


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def parse_length_tsv(tsv_path, tool):
    """Parse a length estimation TSV from telseq or telogator2.

    Returns:
        tuple: (telomere_length_kb, n_measurements)
    """
    df = pd.read_csv(tsv_path, sep="\\t")

    if tool == "telseq":
        if "LENGTH_ESTIMATE" in df.columns:
            numeric_vals = pd.to_numeric(df["LENGTH_ESTIMATE"], errors="coerce").dropna()
            if len(numeric_vals) > 0:
                return numeric_vals.mean(), len(numeric_vals)
        return None, 0

    elif tool == "telogator2":
        if "TL_p75" in df.columns:
            numeric_vals = pd.to_numeric(df["TL_p75"], errors="coerce").dropna()
            if len(numeric_vals) > 0:
                return numeric_vals.median() / 1000.0, len(numeric_vals)
        return None, 0

    return None, 0


def parse_content_tsv(tsv_path):
    """Parse a telomerehunter summary TSV.

    Returns:
        str: telomere content value, or 'NA'
    """
    try:
        df = pd.read_csv(tsv_path, sep="\\t")
        if "tel_content" in df.columns:
            val = pd.to_numeric(df["tel_content"], errors="coerce").dropna()
            if len(val) > 0:
                return str(val.iloc[0])
    except Exception:
        pass
    return "NA"


# Parse length estimation
telomere_length_kb, n_measurements = parse_length_tsv(length_tsv, length_tool)

# Parse content profiling (if available)
content_tool = "NA"
telomere_content = "NA"

if content_tsv and os.path.exists(content_tsv) and os.path.getsize(content_tsv) > 0:
    telomere_content = parse_content_tsv(content_tsv)
    content_tool = "telomerehunter"

# Write summary
summary = pd.DataFrame(
    [
        {
            "sample_id": prefix,
            "length_tool": length_tool,
            "telomere_length_kb": telomere_length_kb if telomere_length_kb is not None else "NA",
            "n_length_measurements": n_measurements,
            "content_tool": content_tool,
            "telomere_content": telomere_content,
        }
    ]
)

summary.to_csv(f"{prefix}_telomere_summary.tsv", sep="\\t", index=False)

# Write versions
versions = {"${task.process}": {"python": platform.python_version(), "pandas": pd.__version__}}
with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
