#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONATOOLS_KTIMPORTTAXONOMY } from '../../../../modules/kronatools/ktimporttaxonomy/main.nf' addParams( options: [:] )

workflow test_kronatools_ktimporttaxonomy {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['generic']['txt']['hello'], checkIfExists: true) ]

    taxonomy = [ file(params.test_data['generic']['txt']['hello'] , checkIfExists: true) ]

    KRONATOOLS_KTIMPORTTAXONOMY ( input, taxonomy )
}
