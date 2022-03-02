#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMETH_INDEX } from '../../../../modules/bwameth/index/main.nf'

workflow test_bwameth_index {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMETH_INDEX ( fasta )
}
