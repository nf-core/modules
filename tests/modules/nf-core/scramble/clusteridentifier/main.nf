#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCRAMBLE_CLUSTERIDENTIFIER } from '../../../../../modules/nf-core/scramble/clusteridentifier/main.nf'

workflow test_scramble_clusteridentifier_bam {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['scramble']['bam'], checkIfExists: true),   
        file(params.test_data['homo_sapiens']['scramble']['bam_bai'], checkIfExists: true),   
        []
    ]

    fasta = []

    SCRAMBLE_CLUSTERIDENTIFIER ( input, fasta )
}

workflow test_scramble_clusteridentifier_cram {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['scramble']['cram'], checkIfExists: true),   
        file(params.test_data['homo_sapiens']['scramble']['cram_crai'], checkIfExists: true),   
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['scramble']['fasta'], checkIfExists: true)

    SCRAMBLE_CLUSTERIDENTIFIER ( input, fasta )
}
