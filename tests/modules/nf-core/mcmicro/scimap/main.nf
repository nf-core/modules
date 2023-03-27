#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MCMICRO_SCIMAP } from '../../../../../modules/nf-core/mcmicro/scimap/main.nf'

workflow test_mcmicro_scimap {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    MCMICRO_SCIMAP ( input )
}
