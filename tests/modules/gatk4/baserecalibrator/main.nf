#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BASERECALIBRATOR } from '../../../../modules/gatk4/baserecalibrator/main.nf' addParams( options: [:] )

workflow test_gatk4_baserecalibrator {
    input     = [ [ id:'test' ], // meta map
                    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
                  ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    sites = file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    sites_tbi = file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)

    GATK4_BASERECALIBRATOR ( input, fasta, fai, dict, [], sites, sites_tbi )
}

workflow test_gatk4_baserecalibrator_intervals {
    input     = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
                ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    intervals = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    sites = file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)
    sites_tbi = file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)

    GATK4_BASERECALIBRATOR ( input, fasta, fai, dict, intervals, sites, sites_tbi )
}

workflow test_gatk4_baserecalibrator_multiple_sites {
    input     = [ [ id:'test' ], // meta map
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
                ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)
    sites =     [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true)
                ]
    sites_tbi = [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true)
                ]

  GATK4_BASERECALIBRATOR ( input, fasta, fai, dict, [], sites, sites_tbi )
}
