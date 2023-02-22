#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTEXT } from '../../../../../modules/nf-core/krona/ktimporttext/main.nf'

workflow test_krona_ktimporttext_multi {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/krona/ktimporttext.txt', checkIfExists: true), // krona default test file
            file(params.test_data['sarscov2']['metagenome']['kraken_report'], checkIfExists: true), //Kraken2 report file
            file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/krona/kaiju_out4krona.txt', checkIfExists: true) // Kaiju output 4 krona
        ]
    ]

    KRONA_KTIMPORTTEXT ( input )
}

workflow test_krona_ktimporttext_single {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file('http://krona.sourceforge.net/examples/text.txt', checkIfExists: true) // krona default test file
        ]
    ]

    KRONA_KTIMPORTTEXT ( input )
}
