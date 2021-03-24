#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_MPILEUP } from '../../../../software/samtools/mpileup/main.nf' addParams( options: [:] )

workflow test_samtools_mpileup {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]
    fasta = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
    
    SAMTOOLS_MPILEUP ( input, fasta )
}
