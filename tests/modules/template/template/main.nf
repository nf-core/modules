#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TEMPLATE_TEMPLATE } from '../../../../modules/template/template/main.nf'

workflow test_template_template {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    TEMPLATE_TEMPLATE ( input )
}
