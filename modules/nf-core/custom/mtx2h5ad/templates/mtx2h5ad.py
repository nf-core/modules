#!/usr/bin/env python

# Nextflow template: select the converter using the module's input_aligner.
import os

# Use a writable mount for the numba cache, including inside containers.
os.environ["NUMBA_CACHE_DIR"] = "."

import glob
import json
import platform

import anndata
import pandas as pd
import scanpy as sc
from anndata import AnnData, concat
from scipy.sparse import csr_matrix


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


def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "pandas": pd.__version__,
            "anndata": anndata.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


def read_cellranger(
    input: str,
):
    adata = sc.read_10x_h5(input)
    adata.var["gene_symbols"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)

    # reorder columns for 10x mtx files
    adata.var = adata.var[["gene_symbols", "feature_types", "genome"]]

    return adata


def read_star(
    input: str,
):
    adata = sc.read_10x_mtx(input)
    adata.layers["count"] = adata.X.copy()

    velocyto_dir = f"velocyto_{input}"
    if os.path.exists(velocyto_dir):
        barcodes = os.path.join(velocyto_dir, "barcodes.tsv.gz")
        features = os.path.join(velocyto_dir, "features.tsv.gz")

        for matrix in ["ambiguous", "spliced", "unspliced"]:
            adata_state = sc.read_mtx(os.path.join(velocyto_dir, f"{matrix}.mtx.gz")).T

            adata_state.obs_names = pd.read_csv(barcodes, header=None, sep="\\t")[0].values
            adata_state.var_names = pd.read_csv(features, header=None, sep="\\t")[0].values

            missing_obs = adata.obs_names[~adata.obs_names.isin(adata_state.obs_names)]
            adata_missing = AnnData(
                X=csr_matrix((len(missing_obs), adata.shape[1])),
                obs=pd.DataFrame(index=missing_obs),
                var=adata_state.var,
            )
            adata_state = concat([adata_state, adata_missing], join="outer")
            adata_state = adata_state[adata.obs_names, adata.var["gene_ids"]].copy()

            adata.layers[matrix] = adata_state.X

    adata.var["gene_symbol"] = adata.var.index
    adata.var.set_index("gene_ids", drop=False, inplace=True)
    return adata


def remove_runtime_metadata(value):
    """Remove measured durations while retaining mapping statistics and provenance."""
    timing_fields = {
        "runtime_seconds",
        "mapping_seconds",
        "busy_nanos",
        "callback_setup_nanos",
        "output_flush_nanos",
    }
    if isinstance(value, dict):
        return {
            key: remove_runtime_metadata(item)
            for key, item in value.items()
            if key not in timing_fields
        }
    if isinstance(value, list):
        return [remove_runtime_metadata(item) for item in value]
    return value


def read_simpleaf(input_data: str):
    # Register the HDF5 filters required by simpleaf output.
    import hdf5plugin

    # open main data
    simpleaf_h5ad_path = f"{input_data}/alevin/quants.h5ad"

    # the simpleaf quant module exports an h5ad file.
    adata = sc.read_h5ad(simpleaf_h5ad_path)
    adata.obs_names = adata.obs["barcodes"].values
    adata.var_names = adata.var["gene_id"].values

    # sort adata column- and row- wise to avoid positional differences
    adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()].copy()

    # Timing measurements and parallel completion order vary between identical runs.
    simpleaf_map_info = remove_runtime_metadata(json.loads(adata.uns["simpleaf_map_info"]))
    adata.uns["simpleaf_map_info"] = json.dumps(simpleaf_map_info, sort_keys=True)
    if "quant_info" in adata.uns:
        quant_info = json.loads(adata.uns["quant_info"])
        if "tiny_cell_resolved_cell_numbers" in quant_info:
            quant_info["tiny_cell_resolved_cell_numbers"] = sorted(quant_info["tiny_cell_resolved_cell_numbers"])
        adata.uns["quant_info"] = json.dumps(quant_info, sort_keys=True)

    return adata


def read_kallisto_matrix(
    matrix: str,
    barcodes: str,
    features: str,
):
    """Load kallisto-formatted mtx files into AnnData."""
    adata = sc.read_mtx(matrix)
    adata.obs_names = pd.read_csv(barcodes, header=None, sep="\\t")[0].values
    adata.var_names = pd.read_csv(features, header=None, sep="\\t")[0].values
    return adata


def combine_kallisto_layers(layers: dict):
    """Align barcode sets, fill missing cells with zeros, and sum counts."""
    all_barcodes = sorted(set().union(*(set(adata.obs_names) for adata in layers.values())))
    aligned = {}
    reference = next(iter(layers.values()))
    for name, adata in layers.items():
        if not adata.var_names.equals(reference.var_names):
            raise ValueError("Kallisto layers must have matching genes in the same order")
        missing = [barcode for barcode in all_barcodes if barcode not in adata.obs_names]
        padding = AnnData(
            X=csr_matrix((len(missing), adata.n_vars)),
            obs=pd.DataFrame(index=missing),
            var=adata.var.copy(),
        )
        aligned[name] = concat([adata, padding], join="outer")[all_barcodes, :].X

    matrices = iter(aligned.values())
    counts = next(matrices).copy()
    for matrix in matrices:
        counts = counts + matrix
    return AnnData(
        X=counts,
        layers=aligned,
        obs=pd.DataFrame(index=all_barcodes),
        var=pd.DataFrame(index=reference.var_names),
    )


def read_kallisto(input_data: str, t2g: str, workflow: str):
    if workflow == "standard":
        adata = read_kallisto_matrix(
            matrix=glob.glob(f"{input_data}/*.mtx")[0],
            barcodes=glob.glob(f"{input_data}/*.barcodes.txt")[0],
            features=glob.glob(f"{input_data}/*.genes.txt")[0],
        )
    elif workflow == "lamanno":
        adata = combine_kallisto_layers({
            state: read_kallisto_matrix(
                matrix=glob.glob(f"{input_data}/{state}*.mtx")[0],
                barcodes=glob.glob(f"{input_data}/{state}*.barcodes.txt")[0],
                features=glob.glob(f"{input_data}/{state}*.genes.txt")[0],
            )
            for state in ("spliced", "unspliced")
        })
    elif workflow == "nac":
        barcodes = glob.glob(f"{input_data}/*.barcodes.txt")[0]
        features = glob.glob(f"{input_data}/*.genes.txt")[0]
        adata = combine_kallisto_layers({
            state: read_kallisto_matrix(
                matrix=glob.glob(f"{input_data}/*{state}.mtx")[0],
                barcodes=barcodes,
                features=features,
            )
            for state in ("nascent", "ambiguous", "mature")
        })
    else:
        raise ValueError(f"Unsupported kallisto workflow: {workflow}")

    txp2gene = pd.read_table(t2g, header=None, names=["gene_id", "gene_symbol"], usecols=[1, 2])
    txp2gene = txp2gene.drop_duplicates(subset="gene_id").set_index("gene_id")
    adata.var = adata.var.join(txp2gene, how="left")
    return adata


def normalize_gene_ids(adata: AnnData):
    """Keep original gene IDs and use unique, version-free IDs as the index."""
    adata.var["gene_versions"] = adata.var.index
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()


def input_to_adata(
    aligner: str,
    input_data: str,
    output: str,
    sample: str,
    input_type: str,
    t2g: str,
    kb_workflow: str,
):
    print(f"Reading {aligner} input for {sample}")
    if aligner in ("cellranger", "cellrangerarc", "cellrangermulti"):
        # cellrangermulti adds a 'sample_' prefix to the matrix filename.
        adata = read_cellranger(glob.glob(f"*{input_type}_feature_bc_matrix.h5")[0])
    elif aligner == "star":
        adata = read_star(input_type)
    elif aligner == "simpleaf":
        adata = read_simpleaf(input_data)
    elif aligner == "kallisto":
        adata = read_kallisto(input_data, t2g, kb_workflow)
    else:
        raise ValueError(f"Unsupported input aligner: {aligner}")

    adata.obs["sample"] = sample
    normalize_gene_ids(adata)
    if aligner == "simpleaf":
        adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()].copy()

    adata.write_h5ad(output)
    print(f"Wrote h5ad file to {output}")
    dump_versions()
    return adata


if __name__ == "__main__":
    input_to_adata(
        aligner="${aligner}",
        input_data="${inputs}",
        output="${prefix}.h5ad",
        sample="${meta.id}",
        input_type="${input_type}",
        t2g="${txp2gene}",
        kb_workflow="${kb_workflow}",
    )
