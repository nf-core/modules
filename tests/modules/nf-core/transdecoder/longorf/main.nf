#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { TRANSDECODER_LONGORF } from "$moduleDir/modules/nf-core/transdecoder/longorf/main.nf"

workflow test_transdecoder_longorf {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    TRANSDECODER_LONGORF ( input )
}
