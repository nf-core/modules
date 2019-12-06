#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../../../nf-core/module_testing/check_process_outputs.nf' params(params)
include '../main.nf' params(params)

// Define input channels
readPaths = [
  ['SRR4238351', ['../../../test-datasets/tools/fastqc/input/SRR4238351_subsamp.fastq.gz']],
  ['SRR4238355', ['../../../test-datasets/tools/fastqc/input/SRR4238355_subsamp.fastq.gz']],
  ['SRR4238359', ['../../../test-datasets/tools/fastqc/input/SRR4238359_subsamp.fastq.gz']],
  ['SRR4238379', ['../../../test-datasets/tools/fastqc/input/SRR4238379_subsamp.fastq.gz']]
]
Channel
  .from(readPaths)
  .map { row -> [ row[0], [ file(row[1][0]) ] ] }
  .set { ch_read_files }

// Run the workflow
workflow {
    fastqc(ch_read_files)
    // .check_output()
}
