#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_REALIGNERTARGETCREATOR } from '../../../../../modules/nf-core/gatk/realignertargetcreator/main.nf'
include { GATK_INDELREALIGNER } from '../../../../../modules/nf-core/gatk/indelrealigner/main.nf'


workflow test_gatk_indelrealigner {

    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    fai = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]
    dict = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    ]

    input_realignertargetcreator = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]

    GATK_REALIGNERTARGETCREATOR ( input_realignertargetcreator, fasta, fai, dict, [[],[]] )

    ch_intervals = GATK_REALIGNERTARGETCREATOR.out.intervals

    ch_bams_indelrealigner = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ])

    ch_input_indelrealigner = ch_bams_indelrealigner.mix(ch_intervals).groupTuple(by: 0).map{ [it[0], it[1][0], it[2], it[1][1] ] }

    GATK_INDELREALIGNER ( ch_input_indelrealigner, fasta, fai, dict, [[],[]] )
}
