#!/usr/bin/env python3

import os
import platform
import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import scanpy.external as sce


# parse_ext_args
# from modules/nf-core/mygene/templates/mygene.py

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
        # input parameters
        parser.add_argument(
            "--columname",
            default="gene_id",
            help="Name of the column where the gene ids are stored in the input file. Default: gene_id",
        )
        # filtering parameters
        parser.add_argument(
            "--species",
            default=None,
            help="Comma separated of common name of the species or taxon ids",
        )
        parser.add_argument(
            "--go_category",
            default=None,
            help="Comma separated list of GO categories to keep. Default: all",
        )
        parser.add_argument(
            "--go_evidence",
            default=None,
            help="Comma separated list of GO evidence codes to keep. Default: all",
        )
        # additional parameters for querymany
        parser.add_argument(
            "--scopes",
            default=None,
            help="Comma separated list of scopes to search for.",
        )
        parser.add_argument(
            "--entrezonly",
            default=False,
            help="When true, the query returns only the hits with valid Entrez gene ids. Default: false.",
        )
        parser.add_argument(
            "--ensemblonly",
            default=False,
            help="When true, the query returns only the hits with valid Ensembl gene ids. Default: False",
        )
        # output parameters
        parser.add_argument(
            "--generate_tsv",
            default=False,
            help="Also generate a tsv file with the gene based information. Default: False",
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

        # check if the arguments are valid
        if args.go_category:
            args.go_category = args.go_category.upper()
            for category in args.go_category.split(","):
                if category not in ["BP", "MF", "CC"]:
                    raise ValueError("The GO category should be one of BP, MF, or CC.")
        if args.go_evidence:
            args.go_evidence = args.go_evidence.upper()

        # Assign args attributes to self attributes
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))

#  ---------------------- new version  ----------------------






param_list = [
    ["hto_data", "${hto_data}"],
    ["priors", "${priors}"],
    ["pre_existing_clusters", "${pre_existing_clusters}"],
    ["number_of_noise_barcodes", "${number_of_noise_barcodes}"],
    ["clustering_data", "${clustering_data}"],
]

# param_df = pd.DataFrame(param_list, columns=["Argument", "Value"])

cell_hashing_data = sc.read_10x_mtx(param_list.hto_data, gex_only=False)
if args.clustering_data is not None:
    trans_data = sc.read_10x_mtx(args.clustering_data)
    trans_data.var_names_make_unique()
    print("--------------------Get data-------------------------------")
    hashsolo.hashsolo(
        cell_hashing_data,
        priors=args.priors,
        clustering_data=trans_data,
        pre_existing_clusters=args.pre_existing_clusters,
        number_of_noise_barcodes=args.number_of_noise_barcodes,
    )
else:
    print("--------------------Get data-------------------------------")
    hashsolo.hashsolo(
        cell_hashing_data,
        priors=args.priors,
        pre_existing_clusters=args.pre_existing_clusters,
        number_of_noise_barcodes=args.number_of_noise_barcodes,
    )
print("--------------------Finished demultiplexing-------------------------------")

print("------------------- Following Files are saved ----------------------------")
print(args.assignmentOutHashSolo + "_res.csv")
print(args.plotOutHashSolo + ".jpg")
print("params.csv")
cell_hashing_data.obs.to_csv(
    args.outputdir + "/" + args.assignmentOutHashSolo + "_res.csv"
)
hashsolo.plot_qc_checks_cell_hashing(cell_hashing_data)
plt.savefig(args.outputdir + "/" + args.plotOutHashSolo + ".jpg", dpi=400)
# param_df.fillna("None", inplace=True)
# param_df.to_csv(args.outputdir + "/params.csv", index=False)



#  ---------------------- old version  ----------------------

adata = sc.read_h5ad("${input_h5ad}")
columns = "${cell_hashing_columns.join(' ')}".split()
columns_str = [str(x) for x in columns]
sce.pp.hashsolo(adata, columns_str, priors=[float(prior) for prior in "${priors.join(',')}".split(',')])

adata.write("${prefix}.h5ad")








# ---------------------- Versions ----------------------

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
