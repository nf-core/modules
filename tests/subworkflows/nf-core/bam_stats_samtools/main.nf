#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_STATS_SAMTOOLS } from '../../../../subworkflows/nf-core/bam_stats_samtools/main'

workflow test_bam_stats_samtools_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_STATS_SAMTOOLS ( input, fasta )
}

workflow test_bam_stats_samtools_paired_end {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BAM_STATS_SAMTOOLS ( input, fasta )
}


workflow test_bam_stats_samtools_paired_end_cram {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BAM_STATS_SAMTOOLS ( input, fasta )
}
