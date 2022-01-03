#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTAXONOMY } from '../../../../modules/krona/ktimporttaxonomy/main.nf'

workflow test_krona_ktimporttaxonomy {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['generic']['txt']['hello'], checkIfExists: true)
    ]
    taxonomy = file(params.test_data['generic']['txt']['hello'], checkIfExists: true)

    KRONA_KTIMPORTTAXONOMY ( input, taxonomy )
}
