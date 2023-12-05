#!/usr/bin/env nextflow

include { MAGUS_ALIGN } from '../../../../../modules/nf-core/magus/align/main.nf'
include { MAGUS_GUIDETREE } from '../../../../../modules/nf-core/magus/guidetree/main.nf'

workflow test_magus_align {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    MAGUS_ALIGN ( input, [ [:], [ ] ] )
}

workflow test_magus_align_with_guide_tree {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    guide_tree = MAGUS_GUIDETREE ( input ).tree

    MAGUS_ALIGN ( input, guide_tree )
}
