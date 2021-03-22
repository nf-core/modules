#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FASTQTOSAM } from '../../../../software/gatk4/fastqtosam/main.nf' addParams( options: [:] )

workflow test_gatk4_fastqtosam_single_end {
    
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ]

    GATK4_FASTQTOSAM ( input )
}
