#!/usr/bin/env nextflow



include { DIAMOND_MAKEDB } from '../../../../modules/diamond/makedb/main.nf'

workflow test_diamond_makedb {

    input = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    DIAMOND_MAKEDB ( input )
}
