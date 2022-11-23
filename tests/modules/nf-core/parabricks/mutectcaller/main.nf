#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_MUTECTCALLER } from '../../../../../modules/nf-core/parabricks/mutectcaller/main.nf'

workflow test_parabricks_mutectcaller {
    
    input = [
        [ tumor_id:'tumour' ],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        [],
        []       
    ]
    normal_input = [
        [],
        [],
        []
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    panel_of_normals = []
    panel_of_normals_index = []

    PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals, panel_of_normals_index )
}

// tn test not passing because of error with readgroups in the test files?
//            [PB Error 2022-Nov-23 15:30:21][src/readers.cpp:1314] ID simulation01 maps to tm_simulation01, but PU is nm_simulation01 (different platform), exiting.                                                                                                                          
// workflow test_parabricks_mutectcaller_tn {
//     input = [
//         [ tumor_id:'tumour' ],
//         file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
//         [],
//         []      
//     ]
//     normal_input = [
//         [ normal_id:'normal' ],
//         file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
//         []
//     ]
//     fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
//     panel_of_normals = []
//     panel_of_normals_index = []

//     PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals, panel_of_normals_index )
// }

workflow test_parabricks_mutectcaller_intervals {
    input = [
        [ tumor_id:'tumour' ],
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
    panel_of_normals_index = []

    PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals, panel_of_normals_index )
}

// PON test not passing because the tool expects all chromosomes in the PON file to 
// be present in the reference genome. 
// workflow test_parabricks_mutectcaller_pon {
//     input = [
//         [ tumor_id:'tumour' ],
//         file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
//         file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
//         file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
//     ]
//     normal_input = [
//         [],
//         [],
//         []
//     ]
//     fasta = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
//     panel_of_normals = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
//     panel_of_normals_index = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

//     PARABRICKS_MUTECTCALLER ( input, normal_input, fasta, panel_of_normals, panel_of_normals_index )
// }
