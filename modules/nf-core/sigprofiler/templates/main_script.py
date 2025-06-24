#!/usr/bin/env python

# Parse arguments

import shlex

def parse_args(x):
    x = x.strip("[]")
    lexer = shlex.shlex(x, posix=True)
    lexer.whitespace = ','
    lexer.whitespace_split = True
    lexer.commenters = ''

    parsed_args = {}
    for item in lexer:
        if ':' not in item:
            continue
        key, val = item.split(':', 1)
        parsed_args[key.strip()] = val.strip().strip("'").strip('"')
    return parsed_args

opt = dict()
opt["prefix"] = "${task.ext.prefix ?: meta.id}"
opt["dataset"] = "${task.ext.dataset ?: meta.dataset}"

opt.update({
    "genome": "",
    "input_type": "matrix",
    "context_type": "96,DINUC,ID",
    "minimum_signatures": 1,
    "maximum_signatures": 25,
    "nmf_replicates": 100,
    "min_nmf_iterations": 10000,
    "max_nmf_iterations": 100000,
    "nmf_test_conv": 10000,
    "seeds": "random",
    "volume": "./",
    "get_all_signature_matrices": True,
    "make_decomposition_plots": True
    })

args_opt = parse_args("${task.ext.args ?: ''}")
for ao_k, ao_v in args_opt.items():
    if ao_k in {"minimum_signatures", "maximum_signatures", "nmf_replicates", "min_nmf_iterations", "max_nmf_iterations", "nmf_test_conv"}:
        opt[ao_k] = int(ao_v)
    elif ao_k in {"get_all_signature_matrices", "make_decomposition_plots"}:
        opt[ao_k] = ao_v.lower() == "true"
    else:
        opt[ao_k] = ao_v


# Script
import os
import shutil
import pandas as pd
import multiprocessing
import subprocess
from importlib.metadata import version


def process_tsv_join(tsv_list):
    patients_tsv = tsv_list.split()
    tables = [pd.read_csv(p, sep='\t', dtype=str) for p in patients_tsv]
    return pd.concat(tables, ignore_index=True)

def input_processing(data, dataset_id, genome):
    new_columns = {'Project': dataset_id, 'Genome': genome, 'Type': "SOMATIC", 'mut_type': "SNP"}
    df = data.assign(**new_columns)
    df['chr'] = df['chr'].astype(str).str[3:]
    df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
    df["ID"] = df["Sample"]
    return df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
    
if __name__ == '__main__':

    dataset_id = opt["dataset"]
    prefix = opt["prefix"]
    genome = opt["genome"]
    input_type = opt["input_type"]
    context_types = opt["context_type"].split(",")
    tsv_list = "${tsv_list.join(' ')}"

    input_path = os.path.join(dataset_id)

    if not os.path.exists(input_path):
            os.mkdir(input_path)
    

    data = process_tsv_join(tsv_list)
    processed = input_processing(data, dataset_id, genome)
    processed.to_csv(f"{input_path}/input_data.txt", sep="\t", index=False)

    # Install genome
    install_genome = f"SigProfilerMatrixGenerator install {genome} -v {opt['volume']}"
    subprocess.run(install_genome, shell=True)
    
    # Mutation counts matrix generation
    generate_count_matrix = (
        f"SigProfilerMatrixGenerator matrix_generator "
        f"-v {opt['volume']} "
        f"{dataset_id} {genome} {input_path}"
    )  
    subprocess.run(generate_count_matrix, shell=True)

    matrix_files = {
            "SBS96": os.path.join("output", "SBS", f"{dataset_id}.SBS96.all"),
            "DBS78": os.path.join("output", "DBS", f"{dataset_id}.DBS78.all"),
            "ID83": os.path.join("output", "ID", f"{dataset_id}.ID83.all")
    }

    context_types = {
        "SBS96": "96",
        "DBS78": "DINUC",
        "ID83": "ID"
    }


    # Run SigProfilerExtractor for each mutation type

    for key, matrix_path in matrix_files.items():
        full_path = os.path.join(dataset_id, matrix_path)

        if not os.path.isfile(full_path):
            continue  # Skip this mutation type

        print(f"Running SigProfilerExtractor on {key} using {full_path}")

        output_dir = f"results/{key}"
        context_type = context_types[key]

        sigprofilerextractor_run = (
            "SigProfilerExtractor sigprofilerextractor "
            "--reference_genome " + genome + " "
            "--context_type " + context_type + " "
            "--minimum_signatures " + str(opt["minimum_signatures"]) + " "
            "--maximum_signatures " + str(opt["maximum_signatures"]) + " "
            "--nmf_replicates " + str(opt["nmf_replicates"]) + " "
            "--seeds " + str(opt["seeds"]) + " "
            "--min_nmf_iterations " + str(opt["min_nmf_iterations"]) + " "
            "--max_nmf_iterations " + str(opt["max_nmf_iterations"]) + " "
            "--nmf_test_conv " + str(opt["nmf_test_conv"]) + " "
            "--get_all_signature_matrices " + str(opt["get_all_signature_matrices"]) + " "
            "--make_decomposition_plots " + str(opt["make_decomposition_plots"]) + " "
            "matrix " + output_dir + " " + full_path
        )  
        subprocess.run(sigprofilerextractor_run, shell=True)
                                                                
                            
    # save the output results
    source_dir = opt["prefix"]

    if not os.path.exists(source_dir):
            os.mkdir(source_dir)

    dest_dir = "results/"
    shutil.copytree(source_dir, "results", dirs_exist_ok=True)

     # Write version
    SigProfilerMatrixGenerator_version = version("SigProfilerMatrixGenerator")
    SigProfilerExtractor_version = version("SigProfilerExtractor")

    with open('versions.yml', 'a') as f:
        f.write('"${task.process}":'+"\\n")
        f.write("    SigProfilerMatrixGenerator: "+SigProfilerMatrixGenerator_version+"\\n")
        f.write("    SigProfilerExtractor: "+SigProfilerExtractor_version+"\\n")

    

