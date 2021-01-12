#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFFREAD } from '../../../software/gffread/main.nf' addParams( options: [:] )

workflow test_gffread {
    def input = []
    input = //[ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/gff/a.gff3", checkIfExists: true) ]
    GFFREAD ( input )
}
