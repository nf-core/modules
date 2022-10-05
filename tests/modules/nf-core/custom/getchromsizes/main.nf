#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUSTOM_GETCHROMSIZES } from '../../../../../modules/nf-core/custom/getchromsizes/main.nf'

workflow test_custom_getchromsizes {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    CUSTOM_GETCHROMSIZES ( input )
}

workflow test_custom_getchromsizes_bgzip {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true) ]

    CUSTOM_GETCHROMSIZES ( input )
}
