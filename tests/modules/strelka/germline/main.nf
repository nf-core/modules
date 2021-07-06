#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRELKA_GERMLINE } from '../../../../modules/strelka/germline/main.nf' addParams( options: [:] )

workflow test_strelka_germline {
    input = [ 
        [ id:'test'], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets = []
    
    STRELKA_GERMLINE ( input, fasta, fai, targets )
}

workflow test_strelka_germline_target_bed {
    input = [ 
        [ id:'test'], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai     = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    targets = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    STRELKA_GERMLINE ( input, fasta, fai, targets )
}

