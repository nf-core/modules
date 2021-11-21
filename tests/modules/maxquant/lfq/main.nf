#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAXQUANT_LFQ } from '../../../../modules/maxquant/lfq/main.nf' addParams( options: [:] )

workflow test_maxquant_lfq {
    
    input = [ [ id:'test' ], // meta map
              file(params.fasta, checkIfExists: true), file(params.paramfile, checkIfExists: true) ]
    rawfiles = file(params.raw)

    MAXQUANT_LFQ ( input, rawfiles.collect())
}
