#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SITEDEPTHTOBAF      } from '../../../../../modules/nf-core/gatk4/sitedepthtobaf/main.nf'
include { GATK4_COLLECTSVEVIDENCE   } from '../../../../../modules/nf-core/gatk4/collectsvevidence/main.nf'

workflow test_gatk4_sitedepthtobaf {

    input = Channel.of([
        [ id:'tumor', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)
    ],
    [
        [ id:'normal', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)
    ]
    )

    vcf = file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true)

    fasta = []
    fasta_fai = []
    dict = []


    GATK4_COLLECTSVEVIDENCE ( input, fasta, fasta_fai, dict )

    sitedepthtobaf_input = GATK4_COLLECTSVEVIDENCE.out.allele_counts.map({ meta, file -> [ [id:'test'], file ]}).groupTuple()

    GATK4_SITEDEPTHTOBAF ( sitedepthtobaf_input, vcf, fasta, dict )
}
