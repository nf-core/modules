#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../../../modules/nf-core/picard/collectmultiplemetrics/main.nf'

workflow test_picard_collectmultiplemetrics {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    fasta = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [[id:'genome'],[]]

    PICARD_COLLECTMULTIPLEMETRICS ( input, fasta, fai )
}

workflow test_picard_collectmultiplemetrics_nofasta {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    PICARD_COLLECTMULTIPLEMETRICS ( input, [[id:'genome'],[]], [[id:'genome'],[]] )
}

workflow test_picard_collectmultiplemetrics_cram {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    fasta = [
        [id:'genome'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id:'genome'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    PICARD_COLLECTMULTIPLEMETRICS ( input, fasta, fai )
}
