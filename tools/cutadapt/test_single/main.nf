#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../main.nf' params(params)

// Define input channels

readPaths = [
  ['SRR4238351', '../../../test-datasets/tools/cutadapt/input/SRR4238351_subsamp.fastq.gz'],
  ['SRR4238355', '../../../test-datasets/tools/cutadapt/input/SRR4238355_subsamp.fastq.gz'],
  ['SRR4238359', '../../../test-datasets/tools/cutadapt/input/SRR4238359_subsamp.fastq.gz'],
  ['SRR4238379', '../../../test-datasets/tools/cutadapt/input/SRR4238379_subsamp.fastq.gz']
]
Channel
  .from(readPaths)
  .map { row -> [ row[0], [ file(row[1]) ] ] }
  .set { ch_read_files }

// Run the workflow
workflow {
    cutadapt(ch_read_files)
}
