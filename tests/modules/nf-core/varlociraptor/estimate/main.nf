#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATE } from '../../../../../modules/nf-core/varlociraptor/estimate/main.nf'

workflow test_varlociraptor_estimate {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_ESTIMATE ( input )
}
