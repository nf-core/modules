#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_SUBSAMPLE_FQ_SALMON } from '../../../../subworkflows/nf-core/fastq_subsample_fq_salmon/main.nf'
include { SALMON_INDEX              } from '../../../../modules/nf-core/salmon/index/main.nf'

//
// Test with paired-end data
//
workflow test_fastq_subsample_fq_salmon {
    input = [
                [ id:'test', single_end:false ], // meta map
                [
                    file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
            ]
    genome_fasta     = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    transcript_fasta = file(params.test_data['homo_sapiens']['genome']['transcriptome_fasta'], checkIfExists: true)
    gtf              = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    SALMON_INDEX ( genome_fasta, transcript_fasta )
    FASTQ_SUBSAMPLE_FQ_SALMON ( input, genome_fasta, transcript_fasta, gtf, SALMON_INDEX.out.index, false )
}
