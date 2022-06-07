#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RHOCALL_ANNOTATE } from '../../../../modules/rhocall/annotate/main.nf'

workflow test_rhocall_annotate {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    RHOCALL_ANNOTATE ( input )
}
