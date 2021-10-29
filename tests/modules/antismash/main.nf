#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANTISMASH } from '../../../modules/antismash/main.nf' addParams( options: [:] )

workflow test_antismash {

    input = [ [ id:'test', single_end:false ], // meta map
              file("https://figshare.com/s/51005a12828ca4390de1", checkIfExists: true) ]

    ANTISMASH ( input )
}
