#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from '../../../../software/samtools/flagstat/main.nf' addParams( options: [:] )

workflow test_samtools_flagstat {
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam.bai", checkIfExists: true) 
            ]
            
    SAMTOOLS_FLAGSTAT ( input )
}
