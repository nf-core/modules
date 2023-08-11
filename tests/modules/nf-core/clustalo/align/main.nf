#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLUSTALO_GUIDETREE } from '../../../../../modules/nf-core/clustalo/guidetree/main.nf'
include { CLUSTALO_ALIGN     } from '../../../../../modules/nf-core/clustalo/align/main.nf'

workflow test_clustalo_align {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]
    
    CLUSTALO_ALIGN ( input, [[:],[]] )
}


workflow test_clustalo_align_with_tree {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    ch_tree = CLUSTALO_GUIDETREE ( input ).tree

    CLUSTALO_ALIGN ( input , ch_tree )

}