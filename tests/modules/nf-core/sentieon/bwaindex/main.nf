#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_BWAINDEX } from '../../../../../modules/nf-core/sentieon/bwaindex/main.nf'

workflow test_sentieon_bwaindex {
    fasta = [ [ id: 'test' ], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    SENTIEON_BWAINDEX ( fasta )
}
