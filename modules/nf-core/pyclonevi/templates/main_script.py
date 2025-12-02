#!/usr/bin/env python

import subprocess

import numpy as np
import pandas as pd


def parse_args(x):
    x = x.replace("[", "")
    x = x.replace("]", "")

    args_list = x.split(", ")
    parsed_args = {i.split(":")[0]: i.split(":")[1] for i in args_list if i.split(":")[1] is not None}
    return parsed_args


opt = dict()
opt["prefix"] = "$task.ext.prefix" if "$task.ext.prefix" == None else "$meta.id"
args_opt = parse_args("$task.ext.args")
for ao_k, ao_v in args_opt.items():
    opt[ao_k] = ao_v

print("$tumour_samples")

# Script #####


# input data preprocessing
def create_pyclone_input(input_data, patient_id, output_data):
    df = pd.read_csv(input_data, sep="\t", header=0)
    # df = df.query('karyotype=="1:1"')
    df["normal_cn"] = 2
    df["patient_id"] = patient_id
    df["mutation_id"] = df.apply(lambda row: f"{patient_id}:{row['chr'][3:]}:{row['from']}:{row['alt']}", axis=1)
    df[["major_cn", "minor_cn"]] = df["karyotype"].str.split(":", expand=True)
    df = df.drop(columns=["karyotype"])
    df = df.rename(columns={"Indiv": "sample_id", "NR": "ref_counts", "NV": "alt_counts", "purity": "tumour_content"})

    column_names = [
        "mutation_id",
        "patient_id",
        "sample_id",
        "ref_counts",
        "alt_counts",
        "normal_cn",
        "major_cn",
        "minor_cn",
        "tumour_content",
        "driver_label",
        "is_driver",
    ]
    column_names_red = list(set(column_names) & set(df.columns))

    df = df.loc[:, column_names_red]
    df.to_csv(output_data, sep="\t", index=False, header=df.columns)


def pyclone_ctree(joint, best_fit, ctree_input):
    ## Read pyclone best fit table
    best_fit_file = best_fit
    df_output = pd.read_csv(best_fit_file, sep="\t")

    ## Read pyclone input table
    joint_table = joint
    df_input = pd.read_csv(joint_table, sep="\t")
    ## Caluclate number of mutations per cluster and add to the subset orginal df_output dataframe

    df_output_small = df_output[["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"]]
    df_output_small["cluster_id"] = "C" + df_output_small["cluster_id"].astype(str)
    df_output_small["nMuts"] = df_output_small.groupby("cluster_id")["mutation_id"].transform("nunique")

    ## Create cluster table
    df_clusters = (
        df_output_small[["sample_id", "cellular_prevalence", "cluster_id"]]
        .drop_duplicates()
        .pivot(index="sample_id", columns="cluster_id", values="cellular_prevalence")
    )
    ## Find the clonal cluster and add the colum 'is.clonal' to identify it

    samples = df_output_small["sample_id"].unique()
    top_clusters = []
    for s in samples:
        max_values = df_clusters.loc[s].max()
        indexes = np.where(df_clusters.loc[s] == max_values)[0]
        top_clusters.append(list(df_clusters.columns[indexes]))
    top_clusters_all = [item for sublist in top_clusters for item in sublist]

    clonal_cluster = max(set(top_clusters_all), key=top_clusters_all.count)

    i = df_output_small[df_output_small["cluster_id"] == clonal_cluster].index

    df_output_small["is.clonal"] = "F"
    df_output_small.loc[i, "is.clonal"] = "T"

    ## Merge the input information about driver genes to the clusters
    if "driver_label" not in df_input.columns or "is_driver" not in df_input.columns:
        df_input["driver_label"] = "NA"
        df_input["is_driver"] = False
        df_input.loc[1, "driver_label"] = ""
        df_input.loc[1, "is_driver"] = True

    df_merged = pd.merge(
        df_output_small,
        df_input[["mutation_id", "sample_id", "patient_id", "driver_label", "is_driver"]],
        on=["mutation_id", "sample_id"],
        how="inner",
    )
    df_final = df_merged.drop_duplicates(subset=["sample_id", "driver_label", "cluster_id", "is.clonal", "is_driver"])
    df_final = df_final.rename(
        columns={
            "cellular_prevalence": "CCF",
            "cluster_id": "cluster",
            "driver_label": "variantID",
            "patient_id": "patientID",
            "is_driver": "is.driver",
        }
    )
    df_final["is.driver"] = df_final["is.driver"].replace({False: "F", True: "T"})

    ## Final ctree input dataframe
    final_jt_file = ctree_input
    df_final.to_csv(final_jt_file, sep="\t", index=False, header=True)


if __name__ == "__main__":
    create_pyclone_input(
        input_data="$rds_join", patient_id="$meta.patient", output_data=opt["prefix"] + "_pyclone_input_all_samples.tsv"
    )

    column_name = "sample_id"
    input_all_samples = pd.read_csv(opt["prefix"] + "_pyclone_input_all_samples.tsv", sep="\t")

    samples_list = "$tumour_samples".replace("[", "")
    samples_list = samples_list.replace("]", "")
    samples_list = samples_list.split(sep=", ")
    input_filtered = input_all_samples[input_all_samples[column_name].isin(samples_list)]
    input_filtered.to_csv(opt["prefix"] + "_pyclone_input.tsv", sep="\t")

    # Run pyclone
    pyclonevi_run = (
        "pyclone-vi fit -i "
        + opt["prefix"]
        + "_pyclone_input.tsv -o "
        + opt["prefix"]
        + "_all_fits.h5 -c "
        + opt["n_cluster"]
    )
    subprocess.run(pyclonevi_run, shell=True)

    pyclonevi_write = (
        "pyclone-vi write-results-file -i " + opt["prefix"] + "_all_fits.h5 -o " + opt["prefix"] + "_best_fit.txt"
    )
    subprocess.run(pyclonevi_write, shell=True)

    # # Pyclone ctree
    pyclone_ctree(
        joint=opt["prefix"] + "_pyclone_input.tsv",
        best_fit=opt["prefix"] + "_best_fit.txt",
        ctree_input=opt["prefix"] + "_cluster_table.csv",
    )

    # Version
    version = (
        subprocess.check_output("pip show pyclone-vi | grep Version | awk '{print \$NF}'", shell=True)
        .decode()
        .split("\\n")[0]
    )
    print(version)
    f = open("versions.yml", "a")
    f.write("$task.process:")
    f.write(f"    pyclone-vi: {version}\\n")
    f.close()
