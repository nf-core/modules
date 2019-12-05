#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../../../nf-core/module_testing/check_process_outputs.nf' params(params)
include '../main.nf' params(params)

// Define input channels
readPaths = [
  ['SRR389222_sub1', ['https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz']],
  ['SRR389222_sub2', ['https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub2.fastq.gz']],
  ['SRR389222_sub3', ['https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub3.fastq.gz']]
]
Channel
  .from(readPaths)
  .map { row -> [ row[0], [row[1][0]]] }
  .set { ch_read_files }

// Run the workflow
workflow {
    fastqc(ch_read_files)
    // .check_output()
}
