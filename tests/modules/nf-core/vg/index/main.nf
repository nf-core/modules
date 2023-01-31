#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VG_INDEX } from '../../../../../modules/nf-core/vg/index/main.nf'

workflow test_vg_index {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VG_INDEX ( input )
}
