#!/usr/bin/env nextflow

include { MAGUS_GUIDETREE } from '../../../../../modules/nf-core/magus/guidetree/main.nf'

workflow test_magus_guidetree {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    MAGUS_GUIDETREE ( input )
}
