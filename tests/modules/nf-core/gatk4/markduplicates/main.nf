#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_MARKDUPLICATES } from '../../../../../modules/nf-core/gatk4/markduplicates/main.nf'
include { GATK4_MARKDUPLICATES as GATK4_MARKDUPLICATES_CRAM } from '../../../../../modules/nf-core/gatk4/markduplicates/main.nf'

workflow test_gatk4_markduplicates {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    GATK4_MARKDUPLICATES ( input, [], [] )
}

workflow test_gatk4_markduplicates_multiple_bams {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
            ] ]

    GATK4_MARKDUPLICATES ( input, [], [] )
}

workflow test_gatk4_markduplicates_multiple_cram_output {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
            ] ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    GATK4_MARKDUPLICATES_CRAM ( input, fasta, fai )
}
