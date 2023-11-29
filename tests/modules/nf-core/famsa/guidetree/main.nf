#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FAMSA_GUIDETREE } from '../../../../../modules/nf-core/famsa/guidetree/main.nf'

workflow test_famsa_guidetree {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    FAMSA_GUIDETREE ( input )
}
