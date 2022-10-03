#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GSTAMA_POLYACLEANUP } from '../../../../modules/gstama/polyacleanup/main.nf'

workflow test_gstama_polyacleanup {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    ]

    GSTAMA_POLYACLEANUP ( input )
}
