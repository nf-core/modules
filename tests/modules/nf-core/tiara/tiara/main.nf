#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TIARA_TIARA } from '../../../../../modules/nf-core/tiara/tiara/main.nf'

workflow test_tiara_tiara {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    TIARA_TIARA ( input )
}
