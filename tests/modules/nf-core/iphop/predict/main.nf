#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IPHOP_DOWNLOAD    } from '../../../../../modules/nf-core/iphop/download/main.nf'
include { IPHOP_PREDICT     } from '../../../../../modules/nf-core/iphop/predict/main.nf'

workflow test_iphop_predict {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)
    ]

    IPHOP_DOWNLOAD ( )

    IPHOP_PREDICT ( input, IPHOP_DOWNLOAD.out.iphop_db )

}
