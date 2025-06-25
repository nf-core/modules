#!/usr/bin/env python3

# Written by Ramprasad Neethiraj and released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.

import argparse
import gzip
import sys
from pathlib import Path
from typing import TextIO


def parse_vep_transcripts(transcripts: list, pli_ind: int) -> list:
    """
    Parse each transcript and return a list of pli values.

    Args:
        transcripts (list): A list of vep transcript annotation
        pli_ind     (int) : Index of pli value in the vep annotation record

    Returns:
        pli_values  (list): list of pli values in the record
    """

    pli_values = []
    for transcript in transcripts:
        vep_fields = transcript.strip().split("|")
        pli_value = vep_fields[pli_ind]
        pli_values.append(pli_value)
    return pli_values


def construct_most_severe_pli_info(line: str, pli_ind: int) -> list:
    """
    Parse gene symbols, find the highest pli value of all gene symbols, add most_severe_pli tag to the info
    field and return a list of modified columns

    Args:
        line        (str) : Vcf record
        pli_ind     (int) : Index of pli value in the vep annotation record

    Returns:
        columns     (list): A list of fields in the vcf record with most severe pli added
                            to the INFO column
    """

    columns = line.strip().split()
    info_fields = columns[7].split(";")
    for field in info_fields:
        if field.startswith("CSQ="):
            transcripts = field.split("CSQ=")[1].split(",")
            break
    pli_values = parse_vep_transcripts(transcripts, pli_ind)
    try:
        pli_max = max(pli_values)
    except ValueError:
        pli_max = ""
    if pli_max:
        columns[7] += ";most_severe_pli={:.2f}".format(float(pli_max))
    return columns


def parse_vep_csq_schema(line: str) -> int:
    """
    Get indices of gene symbol in the annotation

    Args:
        line: INFO line in the vcf header with CSQ information

    Returns:
        pli_ind  (int) : Index of pli value in the vep annotation record
    """
    fields = line.strip().split("Format: ")[1].replace('">', "").split("|")
    pli_ind = fields.index("pLI_gene_value")

    return pli_ind


def write_pli_annotated_vcf(file_in: TextIO, file_out: TextIO):
    """Add most severe pli field to record, and write the record to a vcf file"""
    for line in file_in:
        if line.startswith("#"):
            file_out.write(line)
            if line.startswith("##INFO=<ID=CSQ") and "pLI_gene_value" in line:
                pli_ind = parse_vep_csq_schema(line)
                file_out.write(
                    '##INFO=<ID=most_severe_pli,Number=1,Type=Float,Description="Probability of a gene being loss-of-function intolerant score.">\n'
                )
        else:
            vcf_record = construct_most_severe_pli_info(line, pli_ind)
            file_out.write("\t".join(vcf_record) + "\n")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate vcf with the most severe pli field.",
        epilog="Example: python vcfparser.py --file_in vep.vcf --file_out vep.most_severe_pli.vcf",
    )
    parser.add_argument(
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        help="Vcf file annotated with vep's pli plugin.",
    )
    parser.add_argument(
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Vcf with most_severe_pli annotations added to it.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    if not args.file_in.is_file():
        print(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if (args.file_in.suffix == ".gz") else open
    with open(args.file_out, "w") as out_vcf:
        with opener(args.file_in, "rt") as in_vcf:
            write_pli_annotated_vcf(in_vcf, out_vcf)


if __name__ == "__main__":
    sys.exit(main())
