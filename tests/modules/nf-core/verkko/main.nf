#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VERKKO } from '../../../../modules/nf-core/verkko/main.nf'

workflow test_verkko_pacbio_only {

    input_pacbio = [
        [ id:'test', single_end:true ], // meta map
        [ file("https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz", checkIfExists: true) ]
    ]

     VERKKO ( input_pacbio)
}

