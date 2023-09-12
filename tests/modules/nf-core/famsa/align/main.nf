#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FAMSA_GUIDETREE } from '../../../../../modules/nf-core/famsa/guidetree/main.nf'
include { FAMSA_ALIGN     } from '../../../../../modules/nf-core/famsa/align/main.nf'

workflow test_famsa_align {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]
        
    FAMSA_ALIGN ( input , [[:],[]] )

}


workflow test_famsa_align_with_tree {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    
    ch_tree = FAMSA_GUIDETREE ( input ).tree

    FAMSA_ALIGN ( input , ch_tree )

}
