#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GRAPHMAP2_INDEX } from '../../../../software/graphmap2/index/main.nf' addParams( options: [:] )

workflow test_graphmap2_index {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    GRAPHMAP2_INDEX ( fasta )
}
