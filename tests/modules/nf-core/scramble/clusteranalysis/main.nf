#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCRAMBLE_CLUSTERANALYSIS   } from '../../../../../modules/nf-core/scramble/clusteranalysis/main.nf'
include { SCRAMBLE_CLUSTERIDENTIFIER } from '../../../../../modules/nf-core/scramble/clusteridentifier/main.nf'

workflow test_scramble_clusteranalysis {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['scramble']['bam'], checkIfExists: true),   
        file(params.test_data['homo_sapiens']['scramble']['bam_bai'], checkIfExists: true),   
        []
    ]

    fasta = []
    mei_ref = []

    SCRAMBLE_CLUSTERIDENTIFIER(
        input,
        fasta
    )

    SCRAMBLE_CLUSTERANALYSIS (
        SCRAMBLE_CLUSTERIDENTIFIER.out.clusters,
        fasta,
        mei_ref
    )
}

workflow test_scramble_clusteranalysis_fasta {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['scramble']['cram'], checkIfExists: true),   
        file(params.test_data['homo_sapiens']['scramble']['cram_crai'], checkIfExists: true),   
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['scramble']['fasta'], checkIfExists: true)
    mei_ref = []

    SCRAMBLE_CLUSTERIDENTIFIER(
        input,
        fasta
    )

    SCRAMBLE_CLUSTERANALYSIS (
        SCRAMBLE_CLUSTERIDENTIFIER.out.clusters,
        fasta,
        mei_ref
    )
}