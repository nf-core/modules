#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VT_NORMALIZE } from '../../../../../modules/nf-core/vt/normalize/main.nf'

workflow test_vt_normalize {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
        []
    ]

    fasta = [
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai = [
        [id:'fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VT_NORMALIZE (
        input,
        fasta,
        fai
    )
}

workflow test_vt_normalize_intervals {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    fasta = [
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai = [
        [id:'fai'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VT_NORMALIZE (
        input,
        fasta,
        fai
    )
}

workflow test_vt_normalize_no_fai {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
        []
    ]

    fasta = [
        [id:'fasta'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai = [[],[]]

    VT_NORMALIZE (
        input,
        fasta,
        fai
    )
}
