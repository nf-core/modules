#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_APPLYBQSR } from '../../../../software/gatk4/applybqsr/main.nf' addParams( options: [:] )

workflow test_gatk4_applybqsr {
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/gatk/test.table", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    fai   = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
    dict  = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)

    GATK4_APPLYBQSR ( input, fasta, fai, dict, [] )
}

workflow test_gatk4_applybqsr_intervals {
  input     = [ [ id:'test' ], // meta map
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/gatk/test.table", checkIfExists: true)  
            ]
  fasta     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
  fai       = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
  dict      = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)
  intervals = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)

  GATK4_APPLYBQSR ( input, fasta, fai, dict, intervals )
}
