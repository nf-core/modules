#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SPLITINTERVALS } from '../../../../modules/gatk4/splitintervals/main.nf'

workflow test_gatk4_splitintervals_bed {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    fasta_dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_SPLITINTERVALS ( input, fasta, fasta_fai, fasta_dict)
}

workflow test_gatk4_splitintervals_intervals {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    fasta_dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    GATK4_SPLITINTERVALS ( input, fasta, fasta_fai, fasta_dict)
}