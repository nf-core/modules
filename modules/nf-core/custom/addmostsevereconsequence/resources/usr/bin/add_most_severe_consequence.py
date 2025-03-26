#!/usr/bin/env python3

# Written by Ramprasad Neethiraj and released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.

import argparse
import gzip
import sys
from pathlib import Path
from typing import Tuple, TextIO


def parse_vep_csq_transcripts(
    transcripts: list, allele_ind: int, csq_ind: int, hgnc_ind: int, var_csq: list
) -> Tuple[list, list, list, list]:
    """
    Parse consequences for each transcript and return HGNC IDs, alleles, and their severity rank
    based on the term's ranking in the ensembl consequences list.

    Args:
        transcripts (list): A list of vep transcript annotation
        allele_ind  (int) : Index of the "allele" in the vep annotation record
        csq_ind     (int) : Index of the "Consequence" in the vep annotation record
        hgnc_ind    (int) : Index of the "HGNC_ID" in the vep annotation record
        var_csq     (list): A list of consequence terms ordered by rank

    Returns:
        hgnc_ids     (list): list of hgnc ids in the record
        alleles      (list): list of alleles in the record
        consequences (list): list of consequence terms in the record
        severity     (list): list of consequence term ranks
    """

    consequences = []
    hgnc_ids = []
    severity = []
    alleles = []
    for transcript in transcripts:
        vep_fields = transcript.strip().split("|")
        csq = vep_fields[csq_ind].split("&")[0]
        hgnc_id = vep_fields[hgnc_ind]
        allele = vep_fields[allele_ind].replace("CSQ=", "")
        consequences.append(csq)
        hgnc_ids.append(hgnc_id)
        severity.append(var_csq.index(csq))
        alleles.append(allele)
    return hgnc_ids, alleles, consequences, severity


def construct_most_severe_consequence_info(
    line: str, allele_ind: int, csq_ind: int, hgnc_ind: int, var_csq: list
) -> list:
    """
    Parse consequences for each transcript and return HGNC IDs, alleles, and their severity rank
    based on the term's ranking in the ensembl consequences list.

    Args:
        line        (str) : Vcf record
        allele_ind  (int) : Index of the "allele" in the vep annotation record
        csq_ind     (int) : Index of the "Consequence" in the vep annotation record
        hgnc_ind    (int) : Index of the "HGNC_ID" in the vep annotation record
        var_csq     (list): A list of consequence terms ordered by rank

    Returns:
        columns     (list): A list of fields in the vcf record with most severe consequence added
                            to the INFO column
    """

    columns = line.strip().split()
    info_fields = columns[7].split(";")
    for field in info_fields:
        if field.startswith("CSQ="):
            transcripts = field.split("CSQ=")[1].split(",")
    hgnc_ids, alleles, consequences, severity = parse_vep_csq_transcripts(
        transcripts, allele_ind, csq_ind, hgnc_ind, var_csq
    )
    unique_ids = list(set(hgnc_ids))
    mscsq_anno = []
    for gene_id in unique_ids:
        if gene_id != "":
            indices = find_indices(hgnc_ids, gene_id)
            alleles_sub = [alleles[i] for i in indices]
            consequences_sub = [consequences[i] for i in indices]
            severity_sub = [severity[i] for i in indices]
            most_severe_csq = consequences_sub[severity_sub.index(min(severity_sub))]
            most_severe_allele = alleles_sub[severity_sub.index(min(severity_sub))]
            mscsq_anno.append(gene_id + ":" + most_severe_allele + "|" + most_severe_csq)
    if mscsq_anno:
        columns[7] += ";most_severe_consequence=" + ",".join(mscsq_anno)
    return columns


def find_indices(list_to_check: list, item_to_find: str) -> list:
    """
    Get indices of an element in a list

    Args:
        list_to_check (list)
        item_to_find  (value)

    Returns:
        indices       (list)
    """
    indices = []
    for idx, value in enumerate(list_to_check):
        if value == item_to_find:
            indices.append(idx)
    return indices


def parse_vep_csq_schema(line: str) -> Tuple[int, int, int]:
    """
    Get indices of allele, consequence, and hgnc id in the annotation

    Args:
        line: INFO line in the vcf header with CSQ information

    Returns:
        allele_ind  (int) : Index of the "allele" in the vep annotation record
        csq_ind     (int) : Index of the "Consequence" in the vep annotation record
        hgnc_ind    (int) : Index of the "HGNC_ID" in the vep annotation record
    """
    fields = line.strip().split("Format: ")[1].replace('">', "").split("|")
    allele_ind = fields.index("Allele")
    csq_ind = fields.index("Consequence")
    hgnc_ind = fields.index("HGNC_ID")

    return allele_ind, csq_ind, hgnc_ind


def write_csq_annotated_vcf(file_in: TextIO, file_out: TextIO, var_csq: list):
    """Add most severe consequence field to record, and write the record to a vcf file"""
    for line in file_in:
        if line.startswith("#"):
            file_out.write(line)
            if line.startswith("##INFO=<ID=CSQ"):
                allele_ind, csq_ind, hgnc_ind = parse_vep_csq_schema(line)
                file_out.write(
                    '##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description="Most severe genomic consequence.">\n'
                )
        else:
            mscsq = construct_most_severe_consequence_info(line, allele_ind, csq_ind, hgnc_ind, var_csq)
            file_out.write("\t".join(mscsq) + "\n")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate vcf with the most severe consequence field.",
        epilog="Example: python vcfparser.py --file_in vep.vcf --file_out vep.most_severe_csq.vcf --variant_csq variant_consequence.txt",
    )
    parser.add_argument(
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        help="Vcf file annotated with vep.",
    )
    parser.add_argument(
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Vcf with most_severe_consequence annotations added to it.",
    )
    parser.add_argument(
        "--variant_csq",
        metavar="VARIANT_CSQ",
        type=Path,
        help="Variant consequences ranked by severity",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    if not args.file_in.is_file():
        print(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    if not args.variant_csq.is_file():
        print(f"The given variant consequence file {args.variant_csq} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.variant_csq) as f:
        var_csq = [line.strip() for line in f]
    opener = gzip.open if (args.file_in.suffix == ".gz") else open
    with open(args.file_out, "w") as out_vcf:
        with opener(args.file_in, "rt") as in_vcf:
            write_csq_annotated_vcf(in_vcf, out_vcf, var_csq)


if __name__ == "__main__":
    sys.exit(main())
