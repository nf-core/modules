#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_ADDORREPLACEREADGROUPS as NORMAL_ADDORREPLACEREADGROUPS } from '../../../../../modules/nf-core/picard/addorreplacereadgroups/main.nf'
include { PICARD_ADDORREPLACEREADGROUPS as TUMOUR_ADDORREPLACEREADGROUPS } from '../../../../../modules/nf-core/picard/addorreplacereadgroups/main.nf'
include { SAMTOOLS_INDEX                                                 } from '../../../../../modules/nf-core/samtools/index/main.nf'
include { SENTIEON_TNHAPLOTYPER2                                         } from '../../../../../modules/nf-core/sentieon/tnhaplotyper2/main.nf'

workflow test_tnhaplotyper2 {

    input_normal = [ [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ]
    ]

    input_tumour = [ [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ]
    ]

    NORMAL_ADDORREPLACEREADGROUPS(input_normal)
    TUMOUR_ADDORREPLACEREADGROUPS(input_tumour)

    input = NORMAL_ADDORREPLACEREADGROUPS.out.bam
        .join(NORMAL_ADDORREPLACEREADGROUPS.out.bai)
        .join(TUMOUR_ADDORREPLACEREADGROUPS.out.bam)
        .join(TUMOUR_ADDORREPLACEREADGROUPS.out.bai)
        .map{
            meta, normal_bam, normal_bai, tumour_bam, tumour_bai ->
            [ [id:'test', normal_id:'normal', tumour_id:'tumour'], [normal_bam, tumour_bam], [normal_bai, tumour_bai], [] ]
        }

    dict                  = [ [:], file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true) ]
    fasta                 = [ [:], file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true) ]
    fai                   = [ [:], file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true) ]
    germline_resource     = [ [:], file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true) ]
    germline_resource_tbi = [ [:], file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true) ]
    panel_of_normals      = [ [:], file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true) ]
    panel_of_normals_tbi  = [ [:], file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true) ]

    emit_orientation_data = true
    emit_contamination_data = true

    SENTIEON_TNHAPLOTYPER2(input, dict, fasta, fai, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, emit_orientation_data, emit_contamination_data)

}
