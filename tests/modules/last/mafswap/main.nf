#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LAST_MAFSWAP } from '../../../../modules/last/mafswap/main.nf'

workflow test_last_mafswap {

    input = [ [ id:'contigs.genome' ], // meta map
              file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true) ]

    LAST_MAFSWAP ( input )
}
