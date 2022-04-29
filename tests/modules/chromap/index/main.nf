#!/usr/bin/env nextflow



include { CHROMAP_INDEX } from '../../../../modules/chromap/index/main.nf'

workflow test_chromap_index {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    CHROMAP_INDEX ( input )
}
