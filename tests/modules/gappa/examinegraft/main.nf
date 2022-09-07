#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAPPA_EXAMINEGRAFT } from '../../../../modules/gappa/examinegraft/main.nf'

workflow test_gappa_examinegraft {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GAPPA_EXAMINEGRAFT ( input )
}
