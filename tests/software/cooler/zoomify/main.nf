#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_ZOOMIFY } from '../../../../software/cooler/zoomify/main.nf' addParams( options: ['args':'-r 2,4,8'] )

workflow test_cooler_zoomify {
   input = [ [ id:'test' ], // meta map
            file("https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.asymm.2.cool", checkIfExists: true)]

    COOLER_ZOOMIFY ( input )
}
