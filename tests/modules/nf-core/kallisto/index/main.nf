#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_INDEX } from '../../../../../modules/nf-core/kallisto/index/main.nf'

workflow test_kallisto_index {

    fasta = [
        [ id:'test_fasta' ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    ]

    KALLISTO_INDEX ( fasta )
}
