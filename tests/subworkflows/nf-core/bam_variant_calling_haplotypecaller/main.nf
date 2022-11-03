#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_VARIANT_CALLING_HAPLOTYPECALLER   } from '../../../../subworkflows/nf-core/bam_variant_calling_haplotypecaller/main.nf'
include { GATK4_CALIBRATEDRAGSTRMODEL           } from '../../../../modules/nf-core/gatk4/calibratedragstrmodel/main.nf'

workflow test_bam_variant_calling_haplotypecaller_scatter {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    scatter = Channel.of([
        [id:"test"],
        2
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        scatter,
        fasta,
        fasta_fai,
        dict,
        [],
        [],
        [],
    )
}

workflow test_bam_variant_calling_haplotypecaller_no_scatter {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    scatter = []

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        scatter,
        fasta,
        fasta_fai,
        dict,
        [],
        [],
        [],
    )
}

workflow test_bam_variant_calling_haplotypecaller_scatter_dragstr {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    scatter = Channel.of([
        [id:"test"],
        2
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    strtablefile = file(params.test_data['homo_sapiens']['genome']['genome_strtablefile'], checkIfExists: true)

    GATK4_CALIBRATEDRAGSTRMODEL (
        input,
        fasta,
        fasta_fai,
        dict,
        strtablefile
    )

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        scatter,
        fasta,
        fasta_fai,
        dict,
        GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model,
        [],
        [],
    )
}

workflow test_bam_variant_calling_haplotypecaller_scatter_cram {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    scatter = Channel.of([
        [id:"test"],
        2
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        scatter,
        fasta,
        fasta_fai,
        dict,
        [],
        [],
        [],
    )
}

workflow test_bam_variant_calling_haplotypecaller_scatter_no_scatter {

    input = Channel.of([
        [id:"test"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ],
    [
        [id:"test2"],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ])

    scatter = Channel.of([
        [id:"test"],
        2
    ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    BAM_VARIANT_CALLING_HAPLOTYPECALLER (
        input,
        scatter,
        fasta,
        fasta_fai,
        dict,
        [],
        [],
        [],
    )
}
