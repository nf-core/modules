#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ISLANDPATH } from '../../../../modules/nf-core/islandpath/main.nf'

workflow test_islandpath {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ISLANDPATH ( input )
}
