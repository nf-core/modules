#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_CALLGENERIC } from '../../../../../modules/nf-core/varlociraptor/callgeneric/main.nf'

workflow test_varlociraptor_callgeneric {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_CALLGENERIC ( input )
}
