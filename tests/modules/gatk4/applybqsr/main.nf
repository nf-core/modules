#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_APPLYBQSR } from '../../../../modules/gatk4/applybqsr/main.nf'

workflow test_gatk4_applybqsr {
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
              []
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_APPLYBQSR ( input, fasta, fai, dict )
}

workflow test_gatk4_applybqsr_intervals {
    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
                file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
              ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict  = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

  GATK4_APPLYBQSR ( input, fasta, fai, dict )
}

workflow test_gatk4_applybqsr_cram {
    input = [ [ id:'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
              ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

  GATK4_APPLYBQSR ( input, fasta, fai, dict )
}
