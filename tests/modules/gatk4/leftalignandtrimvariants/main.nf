#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_LEFTALIGNANDTRIMVARIANTS } from '../../../../modules/gatk4/leftalignandtrimvariants/main.nf'

workflow test_gatk4_leftalignandtrimvariants {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict  = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_LEFTALIGNANDTRIMVARIANTS ( input, fasta, fai, dict )
}
