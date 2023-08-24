#!/usr/bin/env python3
from subprocess import run
from pathlib import Path
from textwrap import dedent
import shlex
import re


def chunk_iter(seq, size):
    """iterate over `seq` in chunks of `size`"""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


sample_id = "${meta.id}"

# get fastqs, ordered by path. Files are staged into
#   - "fastq_001/{original_name.fastq.gz}"
#   - "fastq_002/{oritinal_name.fastq.gz}"
#   - ...
# Since we require fastq files in the input channel to be ordered such that a R1/R2 pair
# of files follows each other, ordering will get us a sequence of [R1, R2, R1, R2, ...]
fastqs = sorted(Path(".").glob("fastq_*/*"))
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
    r1.rename(fastq_all / f"{sample_id}_S1_L{i:03d}_R1_001.fastq.gz")
    r2.rename(fastq_all / f"{sample_id}_S1_L{i:03d}_R2_001.fastq.gz")

run(
    # fmt: off
    [
        "cellranger", "count",
        "--id", "${prefix}",
        "--fastqs", str(fastq_all),
        "--transcriptome", "${reference.name}",
        "--localcores", "${task.cpus}",
        "--localmem", "${task.memory.toGiga()}",
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
