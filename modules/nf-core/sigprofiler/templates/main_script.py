#!/usr/bin/env python

import shlex
import os
import shutil
import pandas as pd
import multiprocessing
import subprocess
from importlib.metadata import version

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

# Default options

opt = dict()
opt["prefix"] = "${task.ext.prefix ?: meta.id}"

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
    "make_decomposition_plots": True,
    "download_genome_sigprofiler": True,
    "genome_installed_path": ""
})

# Parse extra args
args_opt = parse_args("${task.ext.args ?: ''}")
for ao_k, ao_v in args_opt.items():
    if ao_k in {"minimum_signatures", "maximum_signatures", "nmf_replicates", "min_nmf_iterations", "max_nmf_iterations", "nmf_test_conv"}:
        opt[ao_k] = int(ao_v)
    elif ao_k in {"get_all_signature_matrices", "make_decomposition_plots", "download_genome_sigprofiler"}:
        val = ao_v.lower()
        if val not in {"true", "false"}:
            raise ValueError(f"Invalid value for '{ao_k}': expected 'true' or 'false', got '{ao_v}'")
        opt[ao_k] = val == "true"
    else:
        opt[ao_k] = ao_v


# Script

def process_tsv_join(tsv_list):
    patients_tsv = tsv_list.split()
    tables = [pd.read_csv(p, sep='\t', dtype=str) for p in patients_tsv]
    return pd.concat(tables, ignore_index=True)

def input_processing(data, prefix, genome):
    new_columns = {'Project': prefix, 'Genome': genome, 'Type': "SOMATIC", 'mut_type': "SNP"}
    df = data.assign(**new_columns)
    df['chr'] = df['chr'].astype(str).str[3:]
    df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
    df["ID"] = df["Sample"]
    return df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
    
if __name__ == '__main__':

    prefix = opt["prefix"]
    genome = opt["genome"]
    input_type = opt["input_type"]
    context_types = opt["context_type"].split(",")
    tsv_list = "${tsv_list.join(' ')}"

    input_path = os.path.join(prefix)
    if not os.path.exists(input_path):
            os.mkdir(input_path)

    data = process_tsv_join(tsv_list)
    processed = input_processing(data, prefix, genome)
    processed.to_csv(f"{input_path}/input_data.txt", sep="\t", index=False)

    # Conditionally install genome or use provided path
    if opt.get("download_genome_sigprofiler", True):
        print(f"Installing genome {genome} via SigProfilerMatrixGenerator...")
        install_genome = f"SigProfilerMatrixGenerator install {genome} -v {opt["volume"]}"
        subprocess.run(install_genome, shell=True)
    else:
        if not opt.get("genome_installed_path"):
            raise ValueError("download_genome_sigprofiler is False but no genome_installed_path was provided.")
        print(f"Using pre-installed genome at: {opt['genome_installed_path']}")
        opt["volume"] 


    # Mutation counts matrix generation
    generate_count_matrix = (
        f"SigProfilerMatrixGenerator matrix_generator "
        f"{prefix} {genome} {input_path} --volume {opt["volume"]}"
    )  
    subprocess.run(generate_count_matrix, shell=True)
    
    matrix_files = {
            "SBS96": os.path.join("output", "SBS", f"{prefix}.SBS96.all"),
            "DBS78": os.path.join("output", "DBS", f"{prefix}.DBS78.all"),
            "ID83": os.path.join("output", "ID", f"{prefix}.ID83.all")
    }

    context_types = {
        "SBS96": "96",
        "DBS78": "DINUC",
        "ID83": "ID"
    }


    # Run SigProfilerExtractor for each mutation type

    for key, matrix_path in matrix_files.items():
        full_path = os.path.join(prefix, matrix_path)

        if not os.path.isfile(full_path):
            continue  # Skip this mutation type

        print(f"Running SigProfilerExtractor on {key} using {full_path}")

        output_dir = f"results/{key}"
        context_type = context_types[key]

        sigprofilerextractor_run = (
            f"SigProfilerExtractor sigprofilerextractor "
            f"--reference_genome  {genome} "
            f"--context_type {context_type} "
            f"--minimum_signatures {opt["minimum_signatures"]} "
            f"--maximum_signatures {opt["maximum_signatures"]} "
            f"--nmf_replicates {opt["nmf_replicates"]} "
            f"--seeds {opt["seeds"]} "
            f"--min_nmf_iterations {opt["min_nmf_iterations"]} "
            f"--max_nmf_iterations {opt["max_nmf_iterations"]} "
            f"--nmf_test_conv {opt["nmf_test_conv"]} "
            f"--get_all_signature_matrices {opt["get_all_signature_matrices"]} "
            f"--make_decomposition_plots {opt["make_decomposition_plots"]} "
            f"matrix {output_dir} {full_path}"
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

    

