#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_MUTECTCALLER } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'

workflow test_parabricks_mutectcaller {
    
    input = [
        [ id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        []
        []       
    ]
    normal_input = [
        [],
        [],
        []
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    panel_of_normals = []
    
    PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals )
}

workflow test_parabricks_mutectcaller_tn {

    input = [
        [ id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        [],
        []       
    ]
    normal_input = [
        [ id:'normal' ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        []
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    panel_of_normals = []

    PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals )
}

workflow test_parabricks_mutectcaller_intervals {
    input = [
        [ id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
    ]
    normal_input = [
        [],
        [],
        []
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    panel_of_normals = []

    PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals )
}

workflow test_parabricks_mutectcaller_pon {
    input = [
        [ id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
    ]
    normal_input = [
        [],
        [],
        []
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)

    PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals )
}
