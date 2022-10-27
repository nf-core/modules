#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { HAPPY_PREPY } from "$moduleDir/modules/nf-core/happy/prepy/main.nf"

workflow test_happy_prepy_vcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome21_indels_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    fasta = Channel.value([
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ])

    HAPPY_PREPY ( input, fasta )
}

workflow test_happy_prepy_gvcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    fasta = Channel.value([
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ])

    HAPPY_PREPY ( input, fasta )
}
