#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_PREPROCESS } from '../../../../../modules/nf-core/varlociraptor/preprocess/main.nf'

workflow test_varlociraptor_preprocess {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_PREPROCESS ( input )
}
