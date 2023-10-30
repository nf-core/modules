#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLUSTALO_GUIDETREE } from '../../../../../modules/nf-core/clustalo/guidetree/main.nf'

workflow test_clustalo_guidetree {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
    ]

    CLUSTALO_GUIDETREE ( input )
}
