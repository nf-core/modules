#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_CLOSEST } from '../../../../../modules/nf-core/bedtools/closest/main.nf'

workflow test_bedtools_closest {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
        [
            file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
        ]
    ]

    fasta_fai = []

    BEDTOOLS_CLOSEST ( input, fasta_fai )
}

workflow test_bedtools_closest_fasta_fai {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ]

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)


    BEDTOOLS_CLOSEST ( input, fasta_fai )
}

workflow test_bedtools_closest_vcf {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'], checkIfExists: true)
    ]

    fasta_fai = []

    BEDTOOLS_CLOSEST ( input, fasta_fai )
}
