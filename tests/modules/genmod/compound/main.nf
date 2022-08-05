#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENMOD_COMPOUND } from '../../../../modules/genmod/compound/main.nf'

workflow test_genmod_compound {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GENMOD_COMPOUND ( input )
}
