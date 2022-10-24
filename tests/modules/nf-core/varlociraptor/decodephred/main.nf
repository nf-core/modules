#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_DECODEPHRED } from '../../../../../modules/nf-core/varlociraptor/decodephred/main.nf'

workflow test_varlociraptor_decodephred {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_bcf'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_DECODEPHRED ( input )
}
