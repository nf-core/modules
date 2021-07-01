#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LAST_DOTPLOT } from '../../../../software/last/dotplot/main.nf' addParams( options: [:] )

workflow test_last_dotplot {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true) ]

    LAST_DOTPLOT ( input, "png" )
}
