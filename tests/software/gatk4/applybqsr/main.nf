#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_APPLYBQSR } from '../../../../software/gatk4/applybqsr/main.nf' addParams( options: [:] )

workflow test_gatk4_applybqsr {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_baserecalibrator_table'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_APPLYBQSR ( input, fasta, fai, dict, [] )
}

workflow test_gatk4_applybqsr_intervals {
    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_baserecalibrator_table'], checkIfExists: true)
              ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    intervals = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

  GATK4_APPLYBQSR ( input, fasta, fai, dict, intervals )
}
