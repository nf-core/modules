#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_CONVERT } from '../../../../modules/bcftools/convert/main.nf'

workflow test_bcftools_convert_gvcf {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
        [],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BCFTOOLS_CONVERT ( input, fasta )
}

workflow test_bcftools_convert_gvcf_bed {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BCFTOOLS_CONVERT ( input, fasta )
}
