#!/home/pollen/anaconda3/bin/python3

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/rnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,file_list,args
    """

    sample_mapping_dict = {}
    with open(file_in, "r", encoding='utf-8-sig') as fin:

        ## Check header
        MIN_COLS = 3
        HEADER = ["sample", "file_list", "args"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )

                num_cols = len([x for x in lspl if x])
                if num_cols < MIN_COLS:
                    print_error(
                        f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                        "Line",
                        line,
                    )

                ## Check sample name entries
                sample, file_list, args = lspl[: len(HEADER)]
                if sample.find(" ") != -1:
                    print(
                        f"WARNING: Spaces have been replaced by underscores for sample: {sample}"
                    )
                    sample = sample.replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)

                ## Check input files
                if file_list:
                    print('file_list: {}'.format(file_list))
                    if file_list.find(" ") != -1:
                        file_list_list = file_list.split(' ')
                    else:
                        file_list_list = [file_list]
                    print(file_list_list)
                    for input_file in file_list_list:
                        print('current input file: {}'.format(input_file))
                        if not os.path.exists(input_file):
                            print_error("Input file does not exist!", "Line", line)
                        if not input_file.endswith(".tiff"):
                            print_error(
                                "Input image file must have extension '.tiff' or '.ome.tiff'!",
                                "Line", line)
                else:
                    print_error("Each run must include at least one input image!", "Line", line)

                ## Check optional arguments
                # TO DO: ashlar will error on bad arguments, so do we need to do any checking here?

                sample_info = [file_list, args]
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                    else:
                        sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(["sample", "file_list", "args"])
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join([f"{sample}"] + val) + "\n")
#                    fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")

    '''
    with open(file_in, "r", encoding='utf-8-sig') as fin:
        with open(file_out, "w") as fout:
            for line in fin:
                fout.write(line)
    '''

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
