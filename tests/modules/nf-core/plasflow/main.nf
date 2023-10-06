#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLASFLOW } from '../../../../modules/nf-core/plasflow/main.nf'

workflow test_plasflow {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    PLASFLOW ( input )
}
