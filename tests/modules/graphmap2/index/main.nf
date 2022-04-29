#!/usr/bin/env nextflow



include { GRAPHMAP2_INDEX } from '../../../../modules/graphmap2/index/main.nf'

workflow test_graphmap2_index {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GRAPHMAP2_INDEX ( fasta )
}
