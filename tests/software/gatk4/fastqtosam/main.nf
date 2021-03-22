#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FASTQTOSAM } from '../../../../software/gatk4/fastqtosam/main.nf' addParams( options: [:] )

workflow test_gatk4_fastqtosam {
    
    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]

    GATK4_FASTQTOSAM ( input )
}
