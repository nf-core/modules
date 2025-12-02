#!/usr/bin/env python3

import argparse
import os
import platform
import shlex

import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import matplotlib
import pandas as pd
import scanpy as sc
import scanpy.external as sce


class Arguments:
    # adopted from mygene module (Suzanne Jin)
    """
    Parses the arguments, including the ones coming from $task.ext.args.
    """

    def __init__(self) -> None:
        self.data = "$data"

        self.use_10x = True
        if self.data.endswith(".h5ad"):
            self.use_10x = False

        cell_hashing_columns = "${cell_hashing_columns.join(' ')}".split()
        self.cell_hashing_columns = [str(x) for x in cell_hashing_columns]

        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

        self.path_assignment = self.prefix + "_assignment_hashsolo.csv"
        # self.path_plot       = self.prefix + "_hashsolo.jpg"
        self.path_h5ad = self.prefix + "_hashsolo.h5ad"
        self.path_params = self.prefix + "_params_hashsolo.csv"
        self.parse_ext_args("$task.ext.args")

    def parse_ext_args(self, args_string: str) -> None:
        """
        It parses the extended arguments.
        """
        # skip when there are no extended arguments
        if args_string == "null":
            args_string = ""

        # Parse the extended arguments
        args_list = shlex.split(args_string)  # Split the string into a list of arguments
        parser = argparse.ArgumentParser()

        parser.add_argument(
            "--priors",
            type=float,
            nargs=3,
            metavar=("NEGATIVE", "SINGLET", "DOUBLET"),
            help="""List of priors for each hypothesis:
            NEGATIVE = prior for negative hypothesis
            SINGLET  = prior for singlet hypothesis
            DOUBLET  = prior for doublet hypothesis""",
            default=[0.01, 0.8, 0.19],
        )

        parser.add_argument(
            "--pre_existing_clusters",
            help="column in cell_hashing_adata.obs for how to break up demultiplexing",
            type=str,
            default=None,
        )
        parser.add_argument(
            "--clustering_data",
            help="input directory containing transcriptomic data in 10x mtx format.",
            type=str,
            default=None,
        )
        parser.add_argument(
            "--number_of_noise_barcodes",
            help="Number of barcodes to use to create noise distribution",
            type=int,
            default=None,
        )

        parser.add_argument(
            "--round_digits",
            help=(
                "Number of decimal places to round numeric values in cell_hashing_data.obs "
                "before saving. If omitted, no rounding is applied."
            ),
            type=int,
            default=None,
        )

        args = parser.parse_args(args_list)

        # Convert "null" values to default values
        # convert true to True and false to False
        for attr in vars(args):
            value = getattr(args, attr)
            if value == "null":
                setattr(args, attr, parser.get_default(attr))
            elif value == "true":
                setattr(args, attr, True)
            elif value == "false":
                setattr(args, attr, False)

        # Assign args attributes to self attributes
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))

    def print_args(self) -> None:
        """
        Print the arguments.
        """
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


if __name__ == "__main__":
    # ------------------------------ parse and print arguments ------------------------------
    args = Arguments()
    args.print_args()

    # ----------------------------------- read input data -----------------------------------
    if args.use_10x:
        cell_hashing_data = sc.read_10x_mtx(args.data, gex_only=False)

        # Move HTO data from .X to .obs columns
        cell_hashing_data.obs[cell_hashing_data.var_names] = pd.DataFrame(
            # Convert sparse matrix to dense array if needed, otherwise use as-is
            cell_hashing_data.X.toarray() if hasattr(cell_hashing_data.X, "toarray") else cell_hashing_data.X,
            index=cell_hashing_data.obs_names,
            columns=cell_hashing_data.var_names,
        )
        args.cell_hashing_columns = list(cell_hashing_data.obs.columns)
    else:
        cell_hashing_data = sc.read_h5ad(args.data)

    if len(args.cell_hashing_columns) == 2:
        # This edge case issue may be fixed in future versions: https://github.com/calico/solo/issues/91
        args.number_of_noise_barcodes = 1

    # -------------------------------------- call hashsolo -----------------------------------
    if args.clustering_data is not None:
        trans_data = sc.read_10x_mtx(args.clustering_data)
        trans_data.var_names_make_unique()
        sce.pp.hashsolo(
            cell_hashing_data,
            cell_hashing_columns=args.cell_hashing_columns,
            priors=args.priors,
            clustering_data=trans_data,
            pre_existing_clusters=args.pre_existing_clusters,
            number_of_noise_barcodes=args.number_of_noise_barcodes,
        )
    else:
        sce.pp.hashsolo(
            cell_hashing_data,
            cell_hashing_columns=args.cell_hashing_columns,
            priors=args.priors,
            pre_existing_clusters=args.pre_existing_clusters,
            number_of_noise_barcodes=args.number_of_noise_barcodes,
        )

    # ------------------------------------- save results -------------------------------------

    # Round numeric values if requested
    if args.round_digits is not None:
        numeric_cols = cell_hashing_data.obs.select_dtypes(include=["number"]).columns
        cell_hashing_data.obs[numeric_cols] = cell_hashing_data.obs[numeric_cols].round(args.round_digits)

    cell_hashing_data.obs.index.name = "Barcode"

    cell_hashing_data.obs.to_csv(args.path_assignment)

    # plotting does not exist for scanpy.external.pp but we couldn't use the solo package
    # sce.pp.hashsolo.plot_qc_checks_cell_hashing(cell_hashing_data)
    # plt.savefig(args.path_plot, dpi=400)

    cell_hashing_data.write(args.path_h5ad)

    param_list = [[key, getattr(args, key)] for key in sorted(vars(args).keys())]

    param_df = pd.DataFrame(param_list, columns=["Argument", "Value"])
    param_df.fillna("None", inplace=True)
    param_df.to_csv(args.path_params, index=False)

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "matplotlib": matplotlib.__version__,
            "pandas": pd.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
