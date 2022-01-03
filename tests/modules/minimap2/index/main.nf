#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIMAP2_INDEX } from '../../../../modules/minimap2/index/main.nf'

workflow test_minimap2_index {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    MINIMAP2_INDEX ( fasta )
}
