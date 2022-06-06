#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_CONSENSUS } from '../../../../modules/bcftools/consensus/main.nf'

workflow test_bcftools_consensus {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BCFTOOLS_CONSENSUS ( input )
}
