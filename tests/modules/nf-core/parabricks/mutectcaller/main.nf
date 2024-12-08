#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_MUTECTCALLER } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'
include { PARABRICKS_MUTECTCALLER as PARABRICKS_MUTECTCALLER_TN } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'
include { PARABRICKS_MUTECTCALLER as PARABRICKS_MUTECTCALLER_TN_INTERVALS } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'
include { PARABRICKS_MUTECTCALLER as PARABRICKS_MUTECTCALLER_PON } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'

workflow test_parabricks_mutectcaller {
    
    input = [
        [ id:'test',  tumor_id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        [], // index only needed if using intervals or tumor-normal calling
        [],
        [],
        []
    ]
    fasta = [
        [ id: 'homo_sapiens' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
        ]
    panel_of_normals = []
    panel_of_normals_index = []

    PARABRICKS_MUTECTCALLER ( input, fasta, panel_of_normals, panel_of_normals_index )
}

// tumor-normal calling not passing because of issue with readgroups in the test files
// [PB Error][src/readers.cpp:1314] ID simulation01 maps to tm_simulation01, but PU is nm_simulation01 (different platform), exiting.
workflow test_parabricks_mutectcaller_tn {
    input = [
        [ id:'test',  tumor_id:'tumour', normal_id:'normal' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        []      
    ]
    fasta = [
         [ id: 'homo_sapiens' ],
         file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
         ]
    panel_of_normals = []
    panel_of_normals_index = []

    PARABRICKS_MUTECTCALLER_TN ( input, fasta, panel_of_normals, panel_of_normals_index )
}

workflow test_parabricks_mutectcaller_intervals {
    input = [
        [ id:'test',  tumor_id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        [],
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
    ]

    fasta = [
        [ id: 'homo_sapiens' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
        ]

    panel_of_normals = []
    panel_of_normals_index = []

    PARABRICKS_MUTECTCALLER ( input, fasta, panel_of_normals, panel_of_normals_index )
}

workflow test_parabricks_mutectcaller_tn_intervals {
    input = [
        [ id:'test',  tumor_id:'tumour', normal_id:'normal'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
    ]
    
    fasta = [
        [ id: 'homo_sapiens' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
        ]
    panel_of_normals = []
    panel_of_normals_index = []

    PARABRICKS_MUTECTCALLER_TN_INTERVALS ( input, fasta, panel_of_normals, panel_of_normals_index )
}

// PON test not passing because the tool expects all chromosomes in the PON file to 
// be present in the reference genome. 
workflow test_parabricks_mutectcaller_pon {
    input = [
        [ id:'test',  tumor_id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        [],
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
    ]

    fasta = [
           [ id: 'homo_sapiens' ],
           file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
           ]

    panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_index = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    PARABRICKS_MUTECTCALLER_PON ( input, fasta, panel_of_normals, panel_of_normals_index )
}
