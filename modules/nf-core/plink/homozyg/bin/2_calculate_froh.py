#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def calculate_froh(hom_file, genome_kb, output_dir):
    # Load ROH .hom file from PLINK
    roh = pd.read_csv(hom_file, sep=r'\s+')

    # Sum ROH length per individual (in kb)
    roh_sum = roh.groupby('FID')['KB'].sum().reset_index()
    roh_sum.columns = ['FID', 'Total_ROH_kb']

    # Compute F_ROH = Total ROH length / Total genome length
    roh_sum['F_ROH'] = roh_sum['Total_ROH_kb'] / genome_kb

    # Save output
    output_file = os.path.join(output_dir, "F_ROH_estimates.csv")
    roh_sum.to_csv(output_file, index=False)
    print(f"[INFO] F_ROH saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate F_ROH from PLINK .hom file")
    parser.add_argument("--roh", required=True, help="Path to PLINK .hom file")
    parser.add_argument("--genome_kb", type=float, default=2500000, help="Total autosomal genome length in kb (default=2.5 Gb)")
    parser.add_argument("--out", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    calculate_froh(args.roh, args.genome_kb, args.out)