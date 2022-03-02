#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LAST_POSTMASK } from '../../../../modules/last/postmask/main.nf'

workflow test_last_postmask {

    input = [ [ id:'contigs.genome' ], // meta map
              file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true) ]

    LAST_POSTMASK ( input )
}
