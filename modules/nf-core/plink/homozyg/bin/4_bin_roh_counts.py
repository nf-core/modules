#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

def bin_roh_counts(hom_path, metadata_path, outdir):
    roh = pd.read_csv(hom_path, sep=r'\s+')
    metadata = pd.read_csv(metadata_path)

    roh['Length_Mb'] = roh['KB'] / 1000
    bins = [1, 2, 4, 8, 16, roh['Length_Mb'].max() + 1]
    labels = ['1-2Mb', '2-4Mb', '4-8Mb', '8-16Mb', '>16Mb']
    roh['ROH_bin'] = pd.cut(roh['Length_Mb'], bins=bins, labels=labels, right=False)

    roh = roh.merge(metadata, on='IID')
    roh_bin_counts = roh.groupby(['IID', 'Breed', 'ROH_bin'], observed=True).size().unstack(fill_value=0)
    roh_bin_counts.to_csv(os.path.join(outdir, "ROH_bins_per_individual.csv"))

    roh_by_breed = roh_bin_counts.groupby('Breed').mean()
    roh_by_breed.to_csv(os.path.join(outdir, "ROH_bins_per_breed.csv"))

    plt.figure(figsize=(10, 6))
    roh_by_breed.T.plot(kind='bar')
    plt.title("Average Number of ROHs per Length Bin by Breed")
    plt.ylabel("Average ROH Count")
    plt.xlabel("ROH Length Bin")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "ROH_Counts_Binned_By_Breed.png"))
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--roh", required=True)
    parser.add_argument("--meta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    bin_roh_counts(args.roh, args.meta, args.out)