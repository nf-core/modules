#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTEXT } from '../../../../modules/krona/ktimporttext/main.nf'

workflow test_krona_ktimporttext {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('http://krona.sourceforge.net/examples/text.txt', checkIfExists: true)
    ]

    KRONA_KTIMPORTTEXT ( input )
}
