#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

sns.set(style="whitegrid")

def visualize_roh(roh_count_file, roh_length_file, outdir, genome_mb=2500):
    os.makedirs(outdir, exist_ok=True)

    roh_counts = pd.read_csv(roh_count_file)
    roh_lengths = pd.read_csv(roh_length_file)

    bins = ['1-2Mb', '2-4Mb', '4-8Mb', '8-16Mb', '>16Mb']

    # A bin with no ROH in the dataset is absent from the step 4/5 tables.
    # Add it as zero so the merge and the per-bin plots still work.
    for b in bins:
        for df in (roh_counts, roh_lengths):
            if b not in df.columns:
                df[b] = 0.0
    # Keep the bins in their usual order, first after IID and Breed.
    roh_lengths = roh_lengths[["IID", "Breed"] + bins + [c for c in roh_lengths.columns if c not in ["IID", "Breed"] + bins]]

    # Merge count + length
    merged = pd.merge(roh_counts, roh_lengths, on=["IID", "Breed"], suffixes=('_count', '_length'))

    # === 1. Correlation plots ===
    print("?? ROH count vs length per bin:")
    scatter_df = pd.DataFrame()
    for bin_label in bins:
        corr = merged[f'{bin_label}_count'].corr(merged[f'{bin_label}_length'])
        print(f"  {bin_label}: r = {corr:.2f}")
        temp = merged[["IID", "Breed", f"{bin_label}_count", f"{bin_label}_length"]].copy()
        temp.columns = ["IID", "Breed", "ROH_Count", "ROH_Length"]
        temp["Bin"] = bin_label
        scatter_df = pd.concat([scatter_df, temp])

    g = sns.lmplot(data=scatter_df, x="ROH_Count", y="ROH_Length", col="Bin", hue="Breed",
                   scatter_kws={"s": 60}, height=4.5, col_wrap=3, ci=None,
                   facet_kws={'sharex': False, 'sharey': False})

    g.set_axis_labels("ROH Count", "ROH Length (Mb)")
    g.fig.subplots_adjust(top=0.88)
    g.fig.suptitle("ROH Count vs Length per Bin")

    for i, bin_label in enumerate(bins):
        subdf = scatter_df[scatter_df["Bin"] == bin_label]
        r, _ = pearsonr(subdf["ROH_Count"], subdf["ROH_Length"])
        g.axes.flatten()[i].text(0.05, 0.95, f"r = {r:.2f}", transform=g.axes.flatten()[i].transAxes, fontsize=10, verticalalignment='top')

    g.savefig(os.path.join(outdir, "ROH_Count_vs_Length_Faceted.png"))
    plt.close()

    # === 2. Boxplot per bin ===
    melted = roh_lengths.melt(id_vars=["IID", "Breed"], var_name="Bin", value_name="ROH_Length")
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=melted, x="Bin", y="ROH_Length", hue="Breed")
    plt.title("ROH Length Distribution per Bin by Breed")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Boxplot_ROH_Length_per_Bin_by_Breed.png"))
    plt.close()

    # === 3. F_ROH per individual ===
    roh_lengths["Total_ROH_Mb"] = roh_lengths[bins].sum(axis=1)
    roh_lengths["F_ROH"] = roh_lengths["Total_ROH_Mb"] / genome_mb
    roh_lengths.to_csv(os.path.join(outdir, "F_ROH_per_individual.csv"), index=False)

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=roh_lengths, x="Breed", y="F_ROH")
    plt.title("F_ROH per Breed")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "F_ROH_by_Breed.png"))
    plt.close()

    # === 4. ROH composition profiles ===
    lengths_norm = roh_lengths.copy()
    lengths_norm[bins] = lengths_norm[bins].div(lengths_norm[bins].sum(axis=1), axis=0)
    lengths_norm_mean = lengths_norm.groupby("Breed")[bins].mean()

    palette = ['#CCB974', '#8172B3', '#4C72B0', '#55A868', '#C44E52']
    lengths_norm_mean.plot(kind='bar', stacked=True, figsize=(10, 6), color=palette)

    plt.ylabel("Proportion of Total ROH Length")
    plt.title("Average ROH Length Composition by Breed")
    plt.xlabel("Breed")
    plt.legend(title="ROH Bin Size", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "ROH_Composition_by_Breed.png"))
    plt.close()

    print("? ROH visualization and statistics complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--roh_counts", required=True)
    parser.add_argument("--roh_lengths", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--genome_mb", type=float, default=2500)
    args = parser.parse_args()
    visualize_roh(args.roh_counts, args.roh_lengths, args.out, args.genome_mb)