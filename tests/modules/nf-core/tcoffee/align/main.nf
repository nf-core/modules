#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FAMSA_GUIDETREE } from '../../../../../modules/nf-core/famsa/guidetree/main.nf'
include { TCOFFEE_ALIGN as TCOFFEE_ALIGN_SEQUENCE } from '../../../../../modules/nf-core/tcoffee/align/main.nf'
include { TCOFFEE_ALIGN as TCOFFEE_ALIGN_WITHTREE } from '../../../../../modules/nf-core/tcoffee/align/main.nf'
 

workflow test_tcoffee_align_sequence {
    
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true)
    ]

    TCOFFEE_ALIGN_SEQUENCE ( input,  [[:],[]],  [[:],[],[]] )
}


workflow test_famsa_align_with_tree {
    
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true)
    ]

    
    ch_tree = FAMSA_GUIDETREE ( input ).tree

    TCOFFEE_ALIGN_WITHTREE ( input , ch_tree,  [[:],[],[]])

}