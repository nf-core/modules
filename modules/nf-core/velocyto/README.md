# Velocyto nf-core module

Velocyto is a library for the analysis of RNA velocity. velocyto.py CLI use
Path(resolve*path=True) and breaks the nextflow logic of symbolic links.
//
Input file must be cellsorted.
To avoid running samtools sort FROM velocyto we need to create a file with
samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam
//
If in the work dir velocyto find a file named EXACTLY cellsorted*[ORIGINALBAMNAME]
it will skip the samtools sort step.
THIS IS WHY I NEED TO HAVE IN work dir a file cellsorted\_ in the same dir
reference: https://velocyto.org/velocyto.py/tutorial/cli.html#notes-on-first-runtime-and-parallelization

NOTE: Optional mask must be passed with ext.args (see tests)
