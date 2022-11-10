#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GERMLINECNVCALLER } from '../../../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'

workflow test_gatk4_germlinecnvcaller {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GATK4_GERMLINECNVCALLER ( input )
}
