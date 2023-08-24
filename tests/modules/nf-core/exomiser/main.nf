#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXOMISER } from '../../../../modules/nf-core/exomiser/main.nf'

workflow test_exomiser {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    EXOMISER ( input )
}
