#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BANDAGE_IMAGE } from '../../../../software/bandage/image/main.nf' addParams( options: [:] )

workflow test_bandage_image {

    def input = []
    input = [ [ id:'B-3106' ], // meta map
              [ file("${launchDir}/tests/data/generic/gfa/B-3106.gfa", checkIfExists: true) ] ]

    BANDAGE_IMAGE ( input )
}
