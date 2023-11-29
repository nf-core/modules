#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DELLY_CALL } from '../../../../../modules/nf-core/delly/call/main.nf'

workflow test_delly_call_bam {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
                    [],
                    [],
                    []
                ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

    DELLY_CALL ( input, fasta, fai )
}

workflow test_delly_call_cram {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                    [],
                    [],
                    []
                ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DELLY_CALL ( input, fasta, fai )
}

workflow test_delly_call_exclude_regions {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                    [],
                    [],
                    file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
                ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DELLY_CALL ( input, fasta, fai )
}

workflow test_delly_call_genotype {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz'], checkIfExists: true),
                    file(params.test_data['homo_sapiens']['illumina']['test_haplotc_cnn_vcf_gz_tbi'], checkIfExists: true),
                    []
                ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DELLY_CALL ( input, fasta, fai )
}
