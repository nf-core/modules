#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GATHERBQSRREPORTS } from '../../../../modules/gatk4/gatherbqsrreports/main.nf'

workflow test_gatk4_gatherbqsrreports {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    GATK4_GATHERBQSRREPORTS ( input )
}
