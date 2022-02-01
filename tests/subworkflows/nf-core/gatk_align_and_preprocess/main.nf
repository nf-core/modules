#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_ALIGN_AND_PREPROCESS } from '../../../../subworkflows/nf-core/gatk_align_and_preprocess/main'

workflow test_gatk_align_and_preprocess_fastq {
    input = [
        [
            [ id:'test', single_end:false, read_group:'"@RG\\tID:Seq01p\\tSM:Seq01\\tPL:ILLUMINA"' ], // meta map
            [
                file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
            ],
            file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
        ]
    ]
    fasta          = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai            = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict           = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    bwaindex       = []
    is_ubam        = false
    sort_order     = "coordinate"
    knownsites     = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    knownsites_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    GATK_ALIGN_AND_PREPROCESS ( input, fasta, fai, dict, bwaindex, is_ubam, sort_order, knownsites, knownsites_tbi )
}

workflow test_gatk_align_and_preprocess_ubam {
    input = [
        [
            [ id:'test', single_end:false, read_group:'"@RG\\tID:Seq01p\\tSM:Seq01\\tPL:ILLUMINA"' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_umi_converted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
        ]
    ]
    fasta          = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai            = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict           = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    bwaindex       = []
    is_ubam        = true
    sort_order     = "coordinate"
    knownsites     = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    knownsites_tbi = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    GATK_ALIGN_AND_PREPROCESS ( input, fasta, fai, dict, bwaindex, is_ubam, sort_order, knownsites, knownsites_tbi )
}
