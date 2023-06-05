#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_REALIGNERTARGETCREATOR } from '../../../../../modules/nf-core/gatk/realignertargetcreator/main.nf'

workflow test_gatk_realignertargetcreator {

    input     = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict  = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK_REALIGNERTARGETCREATOR ( input, fasta, fai, dict, [] )
}
