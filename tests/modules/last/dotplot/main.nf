#!/usr/bin/env nextflow



include { LAST_DOTPLOT } from '../../../../modules/last/dotplot/main.nf'

workflow test_last_dotplot {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true) ]

    LAST_DOTPLOT ( input, "png" )
}
