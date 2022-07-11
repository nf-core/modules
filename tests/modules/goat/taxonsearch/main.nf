#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GOAT_TAXONSEARCH } from '../../../../modules/goat/taxonsearch/main.nf'

workflow test_goat_taxonsearch {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    GOAT_TAXONSEARCH ( input )
}
