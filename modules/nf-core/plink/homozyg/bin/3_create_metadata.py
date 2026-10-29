#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def create_metadata(fam_file, output_dir, name="SheepAnalysis_ROH_metadata.csv"):
    # Read .fam file (6-column format)
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None)
    fam.columns = ["FID", "IID", "Father", "Mother", "Sex", "Phenotype"]

    # Extract metadata: FID = Breed, IID = Sample ID
    metadata = fam[["IID", "FID"]].copy()
    metadata.columns = ["IID", "Breed"]

    # Save metadata
    output_file = os.path.join(output_dir, name)
    metadata.to_csv(output_file, index=False)
    print(f"[INFO] Metadata saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create metadata file from .fam")
    parser.add_argument("--fam", required=True, help="Path to .fam file")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--name", default="SheepAnalysis_ROH_metadata.csv", help="Output file name")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    create_metadata(args.fam, args.out, args.name)