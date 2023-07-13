#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_CONCAT } from '../../../../../modules/nf-core/bcftools/concat/main.nf'

workflow test_bcftools_concat {
    
    input = [ 
        [ id:'test' ], // meta map
        [
          file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true) 
        ],
        [
          file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)
        ]
    ]

    bed = []

    BCFTOOLS_CONCAT ( input, bed )
}

workflow test_bcftools_concat_bed {
    
    input = [ 
        [ id:'test' ], // meta map
        [
          file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true) 
        ],
        [
          file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
          file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true)
        ]
    ]

    bed = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    BCFTOOLS_CONCAT ( input, bed )
}