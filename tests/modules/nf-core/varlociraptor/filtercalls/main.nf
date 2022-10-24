#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_FILTERCALLS } from '../../../../../modules/nf-core/varlociraptor/filtercalls/main.nf'

workflow test_varlociraptor_filtercalls {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_FILTERCALLS ( input )
}
