#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from '../../../../subworkflows/nf-core/gatk_tumor_normal_somatic_variant_calling/main' addParams( [:]  )

workflow test_gatk_tumor_normal_somatic_variant_calling {
    input =         [
                    [ [ id:'test'], // meta map
                    [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true) , file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
                    [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true) , file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
                    ["testN"]
                    ]
                    ]
    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_fasta'],                            checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'],                        checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_dict'],                             checkIfExists: true)
    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'],                    checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'],                checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'],           checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'],       checkIfExists: true)
    interval_file         = file(params.test_data['homo_sapiens']['genome']['genome_interval_list'],                    checkIfExists: true)

    GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING ( input, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, interval_file )
}
