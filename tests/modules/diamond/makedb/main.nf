#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/diamond/makedb/main.nf' addParams( options: [:] )

workflow test_diamond_makedb {

    input = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    DIAMOND_MAKEDB ( input )
}
