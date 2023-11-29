#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_CALIBRATEDRAGSTRMODEL } from '../../../../../modules/nf-core/gatk4/calibratedragstrmodel/main.nf'

workflow test_gatk4_calibratedragstrmodel_bam {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    strtablefile = file(params.test_data['homo_sapiens']['genome']['genome_strtablefile'], checkIfExists: true)

    GATK4_CALIBRATEDRAGSTRMODEL ( input, fasta, fasta_fai, dict, strtablefile )
}

workflow test_gatk4_calibratedragstrmodel_cram {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    strtablefile = file(params.test_data['homo_sapiens']['genome']['genome_strtablefile'], checkIfExists: true)

    GATK4_CALIBRATEDRAGSTRMODEL ( input, fasta, fasta_fai, dict, strtablefile )
}

workflow test_gatk4_calibratedragstrmodel_beds {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    strtablefile = file(params.test_data['homo_sapiens']['genome']['genome_strtablefile'], checkIfExists: true)

    GATK4_CALIBRATEDRAGSTRMODEL ( input, fasta, fasta_fai, dict, strtablefile )
}

workflow test_gatk4_calibratedragstrmodel_gzipped_beds {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    strtablefile = file(params.test_data['homo_sapiens']['genome']['genome_strtablefile'], checkIfExists: true)

    GATK4_CALIBRATEDRAGSTRMODEL ( input, fasta, fasta_fai, dict, strtablefile )
}

