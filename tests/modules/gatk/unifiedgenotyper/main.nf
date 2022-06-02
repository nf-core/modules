#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_UNIFIEDGENOTYPER } from '../../../../modules/gatk/unifiedgenotyper/main.nf'

workflow test_gatk_unifiedgenotyper {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GATK_UNIFIEDGENOTYPER ( input )
}
