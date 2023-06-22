#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VRHYME_VRHYME } from '../../../../../modules/nf-core/vrhyme/vrhyme/main.nf'
include { GUNZIP } from '../../../../../modules/nf-core/gunzip/main.nf'

workflow test_vrhyme_vrhyme {

    bam = [
        [ id:'test', single_end:false ], // meta map
        [
        // file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
        // file(params.test_data['bacteroides_fragilis']['illumina']['test1_2_fastq_gz'], checkIfExists: true),
        // file(params.test_data['bacteroides_fragilis']['illumina']['test2_1_fastq_gz'], checkIfExists: true),
        // file(params.test_data['bacteroides_fragilis']['illumina']['test2_2_fastq_gz'], checkIfExists: true)
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_paired_end_sorted_bam'], checkIfExists: true),
        ]
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        // file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    GUNZIP ( fasta ).gunzip.view()

    VRHYME_VRHYME ( bam, GUNZIP.out.gunzip )
}
