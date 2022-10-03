#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_UNIFIEDGENOTYPER } from '../../../../../modules/nf-core/gatk/unifiedgenotyper/main.nf'

workflow test_gatk_unifiedgenotyper {

    input     = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK_UNIFIEDGENOTYPER ( input, fasta, fai, dict, [], [], [], [])
}
