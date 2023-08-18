#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STITCH } from '../../../../modules/nf-core/stitch/main.nf'

workflow test_stitch {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    STITCH ( input )
}
