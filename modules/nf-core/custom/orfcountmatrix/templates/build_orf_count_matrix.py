#!/usr/bin/env python3
"""Merge per-sample ORF P-site count TSVs into a single ORF x sample matrix.

Input files: any number of per-sample TSVs with three columns:
  sample_id  orf_id  count

The sample_id is identical on every row of a given file (it is prepended
upstream by the per-sample P-site counting step). The column order in the
output matrix is the lexicographic sort of the sample ids encountered.

Output: a single TSV with header `orf_id<tab>sample1<tab>sample2...` and
one row per ORF in the catalogue. ORFs with no counts in a given sample
are filled with 0; ORFs in the catalogue but absent from every sample
appear as a row of zeros. This guarantees the matrix contains every ORF
in the catalogue, not just the subset with non-zero counts in at least
one sample.
"""
import glob
import platform
import sys


def format_yaml_like(data, indent=0):
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def read_orf_ids(path):
    seen = set()
    ordered = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) < 4:
                continue
            orf_id = parts[3]
            if orf_id in seen:
                continue
            seen.add(orf_id)
            ordered.append(orf_id)
    return ordered


def read_counts(path):
    counts = {}
    sample_id = None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) < 3:
                continue
            s, orf, c = parts[0], parts[1], parts[2]
            if sample_id is None:
                sample_id = s
            elif sample_id != s:
                raise ValueError(f"{path}: encountered sample id '{s}' after '{sample_id}'")
            try:
                counts[orf] = counts.get(orf, 0) + int(c)
            except ValueError:
                counts[orf] = counts.get(orf, 0) + int(float(c))
    return sample_id, counts


def build_matrix(count_files, orf_list, output_path):
    orf_ids = read_orf_ids(orf_list)
    per_sample = {}
    for path in count_files:
        sample, counts = read_counts(path)
        if sample is None:
            stem = path.split("/")[-1].split(".")[0]
            sample = stem
            counts = {}
        per_sample[sample] = counts

    samples = sorted(per_sample.keys())
    with open(output_path, "w") as out:
        out.write("orf_id\\t" + "\\t".join(samples) + "\\n")
        for orf in orf_ids:
            row = [orf]
            for s in samples:
                row.append(str(per_sample[s].get(orf, 0)))
            out.write("\\t".join(row) + "\\n")


if __name__ == "__main__":
    count_files = sorted(glob.glob("counts/*"))
    if not count_files:
        sys.stderr.write("No per-sample count files found in counts/\\n")
        sys.exit(1)

    prefix = "$prefix" if "$prefix" != "null" else "orf_psite_counts"
    build_matrix(count_files, "$orf_catalogue_bed12", f"{prefix}.tsv")

    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
