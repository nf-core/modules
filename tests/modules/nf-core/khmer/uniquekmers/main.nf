#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KHMER_UNIQUEKMERS } from '../../../../../modules/nf-core/khmer/uniquekmers/main.nf'

workflow test_khmer_uniquekmers {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    KHMER_UNIQUEKMERS ( input, 50 )
}
