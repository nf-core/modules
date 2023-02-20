#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_MERGE_JASMINE } from '../../../../subworkflows/nf-core/vcf_merge_jasmine/main.nf'

workflow test_vcf_merge_jasmine_default {

    vcfs = Channel.of([
        [id:'test', sample:'sample'],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    ],
    [
        [id:'test2', sample:'sample'],
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
    ])

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    VCF_MERGE_JASMINE (
        vcfs,
        [],
        [],
        fasta,
        fasta_fai,
        [],
        2,
        "sample"
    )
}

workflow test_vcf_merge_jasmine_all_inputs {

    vcfs = Channel.of([
        [id:'test', sample:'sample'],
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    ],
    [
        [id:'test2', sample:'sample'],
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
    ])

    bam = Channel.of([
        [id:'sample', sample:'sample'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ])

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    chr_norm = Channel.of("chr21 21", "chr22 22").collectFile(name:"chr_norm.txt", newLine:true)

    VCF_MERGE_JASMINE (
        vcfs,
        bam,
        [],
        fasta,
        fasta_fai,
        chr_norm,
        2,
        "sample"
    )
}
