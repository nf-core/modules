#!/usr/bin/env python3

"""Run consensus NMF (cNMF; Kotliar et al. 2019) for one rank k.

Input: a features x samples matrix (csv/tsv/txt, optionally gzipped).
Options come from task.ext.args, e.g. '--k 10 --n_iter 200 --seed 42 --prenormalised'.

Two modes:
  default         standard cNMF preprocessing (TPM, high-variance genes, unit-variance
                  scaling), as designed for single-cell counts
  --prenormalised the matrix is used as-is: no TPM, all genes, no scaling, and usages
                  are not refit on TPM (for data that is already normalised)

Outputs (prefix = task.ext.prefix or meta.id):
  <prefix>.feature_loadings.csv    consensus spectra, features x patterns
  <prefix>.sample_factors.csv      usages, samples x patterns
  <prefix>.gene_spectra_score.csv  z-score spectra (program specificity), features x patterns
  <prefix>.stats.csv               k, silhouette (stability), prediction error, settings,
                                   and a summary of the replicate statistics below
  <prefix>.replicate_stats.csv     per NMF replicate: iterations used, whether it hit
                                   --max_nmf_iter, and reconstruction error
  <prefix>.merged_spectra.df.npz   all replicate components before clustering
  <prefix>.clustering.png          cNMF consensus clustergram (if it could be drawn)

Nextflow fills in the dollar-brace placeholders before this script runs, so avoid
literal dollar signs and backslash escapes here.
"""

import argparse
import importlib.metadata
import multiprocessing
import os
import platform
import shlex
import shutil
from concurrent.futures import ProcessPoolExecutor

# These settings must be in place before numpy, scanpy or cNMF are imported, so those
# packages are imported inside the functions below rather than at the top of the file.

# One thread per worker: replicates run in parallel processes instead
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"

# Containers are read-only: numba (used by scanpy) and matplotlib need writable cache folders
for _variable, _folder in (("NUMBA_CACHE_DIR", "numba_cache"), ("MPLCONFIGDIR", "matplotlib_config")):
    os.environ[_variable] = os.path.join(os.getcwd(), _folder)
    os.makedirs(os.environ[_variable], exist_ok=True)

RUN_DIR = "cnmf_run"


def factorize_worker(run_name, worker_index, total_workers):
    """Run one share of the NMF replicates (must be top level for multiprocessing).

    Uses cNMF unchanged, except that each replicate's iterations and reconstruction
    error are recorded: cNMF computes the iteration count but does not keep it.
    """
    import numpy as np
    from cnmf import cNMF
    from sklearn.decomposition import non_negative_factorization

    stats_file = "replicate_stats.worker" + str(worker_index) + ".tsv"

    class RecordingCNMF(cNMF):
        def _nmf(self, data, nmf_kwargs):
            usages, spectra, n_iter = non_negative_factorization(data, **nmf_kwargs)
            dense = data.toarray() if hasattr(data, "toarray") else np.asarray(data)
            error = float(np.linalg.norm(dense - usages @ spectra))
            with open(stats_file, "a") as handle:
                fields = [nmf_kwargs["random_state"], nmf_kwargs["n_components"], n_iter, nmf_kwargs["max_iter"], error]
                handle.write(chr(9).join(str(field) for field in fields) + chr(10))
            return spectra, usages

    RecordingCNMF(output_dir=RUN_DIR, name=run_name).factorize(worker_i=worker_index, total_workers=total_workers)


def collect_replicate_stats(run):
    """Combine the per-worker replicate records and label them with cNMF's replicate numbers."""
    import glob

    import pandas as pd
    from cnmf.cnmf import load_df_from_npz

    columns = ["nmf_seed", "n_components", "n_iter", "max_nmf_iter", "reconstruction_error"]
    records = pd.concat(
        [
            pd.read_csv(path, sep=chr(9), header=None, names=columns)
            for path in sorted(glob.glob("replicate_stats.worker*.tsv"))
        ]
    )
    replicate_params = load_df_from_npz(run.paths["nmf_replicate_parameters"])
    table = replicate_params.merge(records, on=["nmf_seed", "n_components"], how="inner")
    table = table.rename(columns={"iter": "replicate"}).sort_values("replicate")
    table["hit_max_iter"] = table["n_iter"] >= table["max_nmf_iter"]
    return table[["replicate", "nmf_seed", "n_iter", "max_nmf_iter", "hit_max_iter", "reconstruction_error"]]


def read_matrix(file_name):
    """Read a features x samples matrix; feature names in the first column."""
    import pandas as pd

    base_name = file_name[:-3] if file_name.endswith(".gz") else file_name
    separator = "," if base_name.endswith(".csv") else chr(9)
    return pd.read_csv(file_name, sep=separator, index_col=0)


def as_patterns(frame, axis_name):
    """Name the pattern columns Pattern_1..k and the index column, like the cogaps module."""
    frame = frame.copy()
    frame.columns = ["Pattern_" + str(column) for column in frame.columns]
    frame.index.name = axis_name
    return frame


def main():
    import anndata as ad
    import numpy as np
    import pandas as pd
    from cnmf import cNMF
    from cnmf.cnmf import load_df_from_npz

    ## ------------------------------------------------------------------
    ## Options (task.ext.args)
    ## ------------------------------------------------------------------
    parser = argparse.ArgumentParser(prog="cnmf module")
    parser.add_argument("--k", type=int, required=True, help="number of programs")
    parser.add_argument("--n_iter", type=int, default=100, help="number of NMF replicates")
    parser.add_argument("--seed", type=int, default=None, help="random seed")
    parser.add_argument(
        "--density_threshold", type=float, default=0.5, help="outlier filter on local density; 2.0 keeps all components"
    )
    parser.add_argument("--local_neighborhood_size", type=float, default=0.30)
    parser.add_argument("--numgenes", type=int, default=2000, help="high-variance genes (default mode only)")
    parser.add_argument("--max_nmf_iter", type=int, default=1000)
    parser.add_argument("--beta_loss", default="frobenius", choices=["frobenius", "kullback-leibler", "itakura-saito"])
    parser.add_argument(
        "--prenormalised", action="store_true", help="use the matrix as-is: no TPM, all genes, no scaling, no TPM refit"
    )
    ext_args = "${task.ext.args}"
    options = parser.parse_args(shlex.split("" if ext_args == "null" else ext_args))

    if options.k < 2:
        raise SystemExit("--k must be at least 2")
    if options.n_iter < 1:
        raise SystemExit("--n_iter must be at least 1")

    prefix = "${task.ext.prefix}"
    if prefix == "null":
        prefix = "${meta.id}"

    ## ------------------------------------------------------------------
    ## Input: features x samples -> samples x features for cNMF
    ## ------------------------------------------------------------------
    matrix = read_matrix("${matrix}")
    if matrix.isnull().values.any():
        raise SystemExit("The matrix contains missing values")
    if (matrix.values < 0).any():
        raise SystemExit("The matrix contains negative values; NMF needs non-negative input")
    matrix = matrix.loc[matrix.sum(axis=1) > 0]
    if options.k >= matrix.shape[1]:
        raise SystemExit(
            "--k (" + str(options.k) + ") must be smaller than the number of samples (" + str(matrix.shape[1]) + ")"
        )

    samples_by_features = matrix.T
    input_file = "cnmf_input.samples_x_features.tsv"
    samples_by_features.to_csv(input_file, sep=chr(9))
    print(
        "Input:",
        matrix.shape[0],
        "features x",
        matrix.shape[1],
        "samples; k =",
        options.k,
        "; mode =",
        "prenormalised" if options.prenormalised else "standard",
    )

    ## ------------------------------------------------------------------
    ## cNMF: prepare -> factorize (parallel) -> combine -> consensus
    ## ------------------------------------------------------------------
    run = cNMF(output_dir=RUN_DIR, name=prefix)
    if options.prenormalised:
        genes_file = "all_features.txt"
        with open(genes_file, "w") as handle:
            handle.write(chr(10).join(samples_by_features.columns))
        run.prepare(
            counts_fn=input_file,
            components=[options.k],
            n_iter=options.n_iter,
            seed=options.seed,
            tpm_fn=input_file,
            genes_file=genes_file,
            beta_loss=options.beta_loss,
            max_NMF_iter=options.max_nmf_iter,
        )
        # Replace cNMF's unit-variance-scaled matrix with the input, unscaled
        run.save_norm_counts(
            ad.AnnData(
                X=samples_by_features.values.astype(np.float64),
                obs=pd.DataFrame(index=samples_by_features.index),
                var=pd.DataFrame(index=samples_by_features.columns),
            )
        )
    else:
        run.prepare(
            counts_fn=input_file,
            components=[options.k],
            n_iter=options.n_iter,
            seed=options.seed,
            num_highvar_genes=options.numgenes,
            beta_loss=options.beta_loss,
            max_NMF_iter=options.max_nmf_iter,
        )

    total_workers = max(1, min(int("${task.cpus}"), options.n_iter))
    # "spawn" (not the default "fork") so the run shuts down cleanly after the workers finish
    with ProcessPoolExecutor(max_workers=total_workers, mp_context=multiprocessing.get_context("spawn")) as executor:
        list(
            executor.map(
                factorize_worker, [prefix] * total_workers, range(total_workers), [total_workers] * total_workers
            )
        )
    run.combine()
    replicate_stats = collect_replicate_stats(run)

    # Stability and error for this k (the numbers behind cNMF's k-selection plot)
    stats = run.consensus(
        options.k,
        density_threshold=options.density_threshold,
        local_neighborhood_size=options.local_neighborhood_size,
        skip_density_and_return_after_stats=True,
        show_clustering=False,
        close_clustergram_fig=True,
    ).stats
    run.consensus(
        options.k,
        density_threshold=options.density_threshold,
        local_neighborhood_size=options.local_neighborhood_size,
        show_clustering=True,
        close_clustergram_fig=True,
        refit_usage=not options.prenormalised,
    )

    ## ------------------------------------------------------------------
    ## Outputs
    ## ------------------------------------------------------------------
    density_label = str(options.density_threshold).replace(".", "_")
    spectra = load_df_from_npz(run.paths["consensus_spectra"] % (options.k, density_label))
    usages = load_df_from_npz(run.paths["consensus_usages"] % (options.k, density_label))
    scores = load_df_from_npz(run.paths["gene_spectra_score"] % (options.k, density_label))

    as_patterns(spectra.T, "feature").to_csv(prefix + ".feature_loadings.csv")
    as_patterns(usages, "sample").to_csv(prefix + ".sample_factors.csv")
    as_patterns(scores.T, "feature").to_csv(prefix + ".gene_spectra_score.csv")

    # cNMF returns the statistics for one k as a series: k, density threshold, silhouette, prediction error
    stats_row = pd.Series(stats).to_dict() if not isinstance(stats, pd.DataFrame) else stats["stats"].to_dict()
    replicate_stats.to_csv(prefix + ".replicate_stats.csv", index=False)
    errors = replicate_stats["reconstruction_error"]
    pd.DataFrame(
        [
            {
                "k": options.k,
                "n_iter": options.n_iter,
                "density_threshold": options.density_threshold,
                "local_neighborhood_size": options.local_neighborhood_size,
                "prenormalised": options.prenormalised,
                "silhouette": stats_row["silhouette"],
                "prediction_error": stats_row["prediction_error"],
                "max_nmf_iter": options.max_nmf_iter,
                "mean_nmf_iterations": replicate_stats["n_iter"].mean(),
                "max_nmf_iterations": replicate_stats["n_iter"].max(),
                "replicates_hit_max_iter": int(replicate_stats["hit_max_iter"].sum()),
                "reconstruction_error_min": errors.min(),
                "reconstruction_error_max": errors.max(),
                "reconstruction_error_std": errors.std(ddof=0),
                "reconstruction_error_cv": errors.std(ddof=0) / errors.mean() if errors.mean() else float("nan"),
            }
        ]
    ).to_csv(prefix + ".stats.csv", index=False)

    shutil.copy(run.paths["merged_spectra"] % options.k, prefix + ".merged_spectra.df.npz")
    clustering_plot = run.paths["clustering_plot"] % (options.k, density_label)
    if os.path.exists(clustering_plot):
        shutil.copy(clustering_plot, prefix + ".clustering.png")

    print(
        "Replicates hitting --max_nmf_iter:",
        int(replicate_stats["hit_max_iter"].sum()),
        "of",
        len(replicate_stats),
    )
    print(
        "Silhouette:",
        round(float(stats_row["silhouette"]), 4),
        "| prediction error:",
        round(float(stats_row["prediction_error"]), 2),
    )

    ## versions.yml
    with open("versions.yml", "w") as handle:
        handle.write('"${task.process}":' + chr(10))
        for package in ("cnmf", "scikit-learn", "numpy", "pandas"):
            handle.write("    " + package + ": " + importlib.metadata.version(package) + chr(10))
        handle.write("    python: " + platform.python_version() + chr(10))


if __name__ == "__main__":
    main()
