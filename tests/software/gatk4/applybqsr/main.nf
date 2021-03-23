#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_APPLYBQSR } from '../../../../software/gatk4/applybqsr/main.nf' addParams( options: [:] )

workflow test_gatk4_applybqsr {

  def input = []
  input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
            file("${launchDir}/tests/data/genomics/sarscov2/table/test.table", checkIfExists: true) ]
  fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
  fastaidx = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
  dict = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.dict", checkIfExists: true)

    GATK4_APPLYBQSR ( input, fasta, fastaidx, dict, [] )
}

workflow test_gatk4_applybqsr_intervals {

  def input = []
  input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true),
            file("${launchDir}/tests/data/genomics/sarscov2/table/test.table", checkIfExists: true)  ]
  fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
  fastaidx = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
  dict = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.dict", checkIfExists: true)
  intervals = file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true)

  GATK4_APPLYBQSR ( input, fasta, fastaidx, dict, intervals)
}
