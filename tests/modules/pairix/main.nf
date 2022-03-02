#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRIX } from '../../../modules/pairix/main.nf'

workflow test_pairix {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://raw.githubusercontent.com/4dn-dcic/pairix/master/samples/4dn.bsorted.chr21_22_only.nontriangle.pairs.gz", checkIfExists: true) ]

    PAIRIX ( input )
}
