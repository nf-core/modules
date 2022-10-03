#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES_SPARK } from '../../../../../modules/nf-core/gatk4/markduplicatesspark/main.nf'
include { GATK4_MARKDUPLICATES_SPARK as GATK4_MARKDUPLICATES_SPARK_CRAM } from '../../../../../modules/nf-core/gatk4/markduplicatesspark/main.nf'
include { GATK4_MARKDUPLICATES_SPARK as GATK4_MARKDUPLICATES_SPARK_METRICS } from '../../../../../modules/nf-core/gatk4/markduplicatesspark/main.nf'

workflow test_gatk4_markduplicates_spark {
    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MARKDUPLICATES_SPARK ( input, fasta, fai, dict )
}

// chr 22
workflow test_gatk4_markduplicates_spark_multiple_bams {
    input = [ [ id:'test', single_end:false ], // meta map
            [   file(params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true)
            ] ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MARKDUPLICATES_SPARK ( input, fasta, fai, dict )
}

// chr 22
workflow test_gatk4_markduplicates_spark_multiple_bams_cram_out {
    input = [ [ id:'test', single_end:false ], // meta map
            [   file(params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true)
            ] ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MARKDUPLICATES_SPARK_CRAM ( input, fasta, fai, dict )
}

// chr 22
workflow test_gatk4_markduplicates_spark_multiple_bams_metrics {
    input = [ [ id:'test', single_end:false ], // meta map
            [   file(params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true)
            ] ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_MARKDUPLICATES_SPARK_METRICS ( input, fasta, fai, dict )
}
