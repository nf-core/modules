#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_FILTERVARIANTTRANCHES } from '../../../../modules/gatk4/filtervarianttranches/main.nf'

workflow test_gatk4_filtervarianttranches {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GATK4_FILTERVARIANTTRANCHES ( input )
}
