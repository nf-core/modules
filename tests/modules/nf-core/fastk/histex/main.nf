#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { FASTK_FASTK } from "$moduleDir/modules/nf-core/fastk/fastk/main.nf"
include { FASTK_HISTEX } from "$moduleDir/modules/nf-core/fastk/histex/main.nf"

workflow test_fastk_histex {

    input = [
        [ id:'test' , single_end: true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
    FASTK_HISTEX ( FASTK_FASTK.out.hist )
}
