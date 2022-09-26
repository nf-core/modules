#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../../modules/nf-core/picard/collectmultiplemetrics/main.nf'

workflow test_picard_collectmultiplemetrics {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    PICARD_COLLECTMULTIPLEMETRICS ( input, fasta, [] )
}

workflow test_picard_collectmultiplemetrics_nofasta {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    PICARD_COLLECTMULTIPLEMETRICS ( input, [], [] )
}

workflow test_picard_collectmultiplemetrics_cram {
    input = [
                [ id:'test', single_end:false ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true)
            ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    PICARD_COLLECTMULTIPLEMETRICS ( input, fasta, fai )
}
