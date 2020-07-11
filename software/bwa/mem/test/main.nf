#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../../../../nf-core/module_testing/check_process_outputs.nf' params(params)
include '../main.nf' params(params)

reads = '../../../../test-datasets/tools/bwa/mem/reads/*_R{1,2}_001.fastq.gz'
index = '../../../../test-datasets/tools/bwa/mem/index/H3N2.{amb,ann,bwt,pac,sa}'
prefix = 'H3N2'

workflow {
  read_input=Channel.fromFilePairs(reads)
  bwa_mem(read_input,file(index),prefix)
}
