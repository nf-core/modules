#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../../modules/nf-core/diamond/makedb/main.nf'

workflow test_diamond_makedb {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    DIAMOND_MAKEDB ( input )
}
