#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMCOPY_GENERATEMAP } from '../../../../modules/hmmcopy/generatemap/main.nf'

workflow test_hmmcopy_generatemap {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    HMMCOPY_GENERATEMAP ( fasta )
}
