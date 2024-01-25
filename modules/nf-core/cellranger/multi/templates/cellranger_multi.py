#!/usr/bin/env python3
"""
Automatically rename staged files for input into cellranger multi and run it.

Copyright (c) Felipe Almeida 2024 - MIT License
"""
from subprocess import run
from pathlib import Path
from textwrap import dedent
import shlex
import re


def chunk_iter(seq, size):
    """iterate over `seq` in chunks of `size`"""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


work_dir  = str(Path.cwd())
sample_id = "${prefix}"

# get fastqs, ordered by path. Files are staged into
#   - "fastq_001/{original_name.fastq.gz}"
#   - "fastq_002/{original_name.fastq.gz}"
#   - ...
# Since we require fastq files in the input channel to be ordered such that a R1/R2 pair
# of files follows each other, ordering will get us a sequence of [R1, R2, R1, R2, ...]
fastqs = sorted(Path(".").glob("fastqs/*/[!EMPTY]*"))
assert len(fastqs) % 2 == 0

# target directory in which the renamed fastqs will be placed
fastq_all = Path("./fastq_all")
fastq_all.mkdir(exist_ok=True)

# Match R1 in the filename, but only if it is followed by a non-digit or non-character
# match "file_R1.fastq.gz", "file.R1_000.fastq.gz", etc. but
# do not match "SRR12345", "file_INFIXR12", etc
filename_pattern =  r'([^a-zA-Z0-9])R1([^a-zA-Z0-9])'

for i, (r1, r2) in enumerate(chunk_iter(fastqs, 2)):
    # double escapes are required because nextflow processes this python 'template'
    if re.sub(filename_pattern, r'\\1R2\\2', r1.name) != r2.name:
        raise AssertionError(
            dedent(
                f"""\
                We expect R1 and R2 of the same sample to have the same filename except for R1/R2.
                This has been checked by replacing "R1" with "R2" in the first filename and comparing it to the second filename.
                If you believe this check shouldn't have failed on your filenames, please report an issue on GitHub!

                Files involved:
                    - {r1}
                    - {r2}
                """
            )
        )

    # get sub-dir of fastqs
    subdir = str(r1).split('/')[-2]

    # target directory in which the renamed fastqs will be placed
    final_dir = Path(f"./fastq_all/{subdir}")
    final_dir.mkdir(exist_ok=True)

    # rename
    r1.rename(fastq_all / subdir / f"{sample_id}_S1_L{i:03d}_R1_001.fastq.gz")
    r2.rename(fastq_all / subdir / f"{sample_id}_S1_L{i:03d}_R2_001.fastq.gz")

#
# fix relative paths from main.nf
#
work_dir = str(Path.cwd()) + '/'
gex_reference_path = "${gex_reference_path}".replace('./', work_dir)
frna_probeset = "${frna_probeset}".replace('./', work_dir)
target_panel = "${target_panel}".replace('./', work_dir)
cmo_reference_path = "${cmo_reference_path}".replace('./', work_dir)
cmo_barcode_path = "${cmo_barcode_path}".replace('./', work_dir)
fb_reference_path = "${fb_reference_path}".replace('./', work_dir)
vdj_reference_path = "${vdj_reference_path}".replace('./', work_dir)
primer_index = "${primer_index}".replace('./', work_dir)
fastq_gex = "${fastq_gex}".replace('./', work_dir)
fastq_vdj = "${fastq_vdj}".replace('./', work_dir)
fastq_antibody = "${fastq_antibody}".replace('./', work_dir)
fastq_beam = "${fastq_beam}".replace('./', work_dir)
fastq_crispr = "${fastq_crispr}".replace('./', work_dir)
fastq_cmo = "${fastq_cmo}".replace('./', work_dir)

#
# generate config file for cellranger multi
#
config_txt = f"""${include_gex}
{gex_reference_path}
{frna_probeset}
${gex_options_filter_probes}
${gex_options_r1_length}
${gex_options_r2_length}
${gex_options_chemistry}
${gex_options_expect_cells}
${gex_options_force_cells}
${gex_options_no_secondary}
${gex_options_no_bam}
${gex_options_check_library_compatibility}
{target_panel}
${gex_options_no_target_umi_filter}
${gex_options_include_introns}
${cmo_options_min_assignment_confidence}
{cmo_reference_path}
{cmo_barcode_path}

${include_fb}
{fb_reference_path}
${fb_options_r1_length}
${fb_options_r2_length}

${include_vdj}
{vdj_reference_path}
{primer_index}
${vdj_options_r1_length}
${vdj_options_r2_length}

[libraries]
fastq_id,fastqs,lanes,feature_types
{fastq_gex}
{fastq_vdj}
{fastq_antibody}
{fastq_beam}
{fastq_crispr}
{fastq_cmo}
"""

with open("${config}", 'w') as file:
    file.write(config_txt)

#
# check the extra data that is included
#
if len("${include_cmo}") > 0:
    with open("${config}", 'a') as file, open("${cmo_barcodes}", 'r') as input_conf:
        file.write("${include_cmo}\\n")
        for line in input_conf:
            file.write(line + "\\n")

if len("${include_beam}") > 0:
    with open("${config}", 'a') as file, open("${beam_csv_text}", 'r') as input_conf, open("${beam_antigen_csv}", 'r') as input_csv:
        file.write("${include_beam}\\n")
        for line in input_conf:
            file.write(line + "\\n")
        file.write("[feature]\\n")
        for line in input_csv:
            file.write(line + "\\n")

if len("${include_frna}") > 0:
    with open("${config}", 'a') as file, open("${frna_csv_text}", 'r') as input_conf:
        file.write("${include_frna}\\n")
        for line in input_conf:
            file.write(line + "\\n")

# Remove blank lines from config file
with open('tmp.txt', 'w') as file:
    file.write( run(['grep', '-v', '-e', '^[[:space:]]*\$', "${config}"], capture_output=True, text=True).stdout )
run(['mv', 'tmp.txt', "${config}"])

#
# run cellranger multi
#
run(
    # fmt: off
    [
        'cellranger', 'multi',
        f"--id={sample_id}",
        "--csv=${config}",
        "--localcores=${task.cpus}",
        "--localmem=${task.memory.toGiga()}",
        *shlex.split("""${args}""")
    ],
    # fmt: on
    check=True,
)

# Output version information
version = run(
    ["cellranger", "-V"],
    text=True,
    check=True,
    capture_output=True,
).stdout.replace("cellranger cellranger-", "")

# alas, no `pyyaml` pre-installed in the cellranger container
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'  cellranger: "{version}"\\n')
