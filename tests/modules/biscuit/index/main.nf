#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_INDEX } from '../../../../modules/biscuit/index/main.nf'

workflow test_biscuit_index {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX ( fasta )
}
