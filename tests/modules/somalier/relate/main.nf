#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_RELATE } from '../../../../modules/somalier/relate/main.nf'

input_sample_ch = Channel
.fromPath( params.inputCsv )
.splitCsv( header:true )
.map { row -> tuple( [id: row.patient], file(row.somalierextract) ) }

workflow test_somalier_relate {

    /*input = [
        [ id:'all', single_end:false ], // meta map
        [file("/Users/asma/Downloads/normal.somalier", checkIfExists: true),
        file("/Users/asma/Downloads/tumour.somalier", checkIfExists: true)]
    ]
    */

    SOMALIER_RELATE (input_sample_ch.collect { it[1] })
}
