#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_APPLYBQSR } from '../../../../../modules/nf-core/parabricks/applybqsr/main.nf'

workflow test_parabricks_applybqsr {
    
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        [], // index not needed unless using intervals
        file(params.test_data['sarscov2']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
        []
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]

    PARABRICKS_APPLYBQSR ( input, fasta )
}

workflow test_parabricks_applybqsr_intervals {
    
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_baserecalibrator_table'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]
        
    PARABRICKS_APPLYBQSR ( input, fasta )
}
