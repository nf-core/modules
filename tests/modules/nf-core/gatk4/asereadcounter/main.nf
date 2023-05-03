#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_ASEREADCOUNTER } from '../../../../../modules/nf-core/gatk4/asereadcounter/main.nf'

workflow test_gatk4_asereadcounter {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]
    vcf = [ [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    intervals = file(params.test_data['sarscov2']['genome']['targets_interval_list'], checkIfExists: true)

    GATK4_ASEREADCOUNTER ( input, vcf, fasta, fai, dict, intervals )
}
