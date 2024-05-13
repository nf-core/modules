#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VG_CONSTRUCT } from '../../../../../modules/nf-core/vg/construct/main.nf'

workflow test_vg_construct_vcf {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        []
    ]

    fasta = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]


    VG_CONSTRUCT (
        input,
        fasta,
        fasta_fai
    )
}

workflow test_vg_construct_msa {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa'], checkIfExists: true),
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    fasta = [[],[]]
    fasta_fai = [[],[]]

    VG_CONSTRUCT (
        input,
        fasta,
        fasta_fai
    )
}

workflow test_vg_construct_insertions_fasta_vcf {

    input = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['pangenome']['pangenome_fa'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'test'], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    ]

    VG_CONSTRUCT (
        input,
        fasta,
        fasta_fai
    )
}
