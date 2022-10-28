#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_VARIANT_CALLING_HAPLOTYPECALLER } from '../../../../subworkflows/nf-core/bam_variant_calling_haplotypecaller/main.nf'

workflow test_bam_variant_calling_haplotypecaller_scatter {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        fasta,
        fasta_fai,
        dict,
        [],
        [],
        [],
        "scatter",
        2
    )
}

workflow test_bam_variant_calling_haplotypecaller_no_scatter {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        fasta,
        fasta_fai,
        dict,
        [],
        [],
        [],
        "no_scatter",
        []
    )
}
