#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BACPHLIP } from '../../../../modules/nf-core/bacphlip/main.nf'

workflow test_bacphlip {

    input = [
        [ id:'test', single_end:false ], // meta map
        file('https://bitbucket.org/srouxjgi/iphop/raw/d27b6bbdcd39a6a1cb8407c44ccbcc800d2b4f78/test/test_input_phages.fna', checkIfExists: true)
    ]

    BACPHLIP ( input )
}
