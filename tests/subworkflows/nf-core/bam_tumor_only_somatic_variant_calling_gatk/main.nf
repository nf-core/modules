#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK } from '../../../../subworkflows/nf-core/bam_tumor_only_somatic_variant_calling_gatk/main'

workflow test_bam_tumor_only_somatic_variant_calling_gatk {
    input = Channel.of([
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'],     checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        []
    ])

    fasta                 = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'],                      checkIfExists: true))
    fai                   = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'],                  checkIfExists: true))
    dict                  = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_dict'],                       checkIfExists: true))
    germline_resource     = Channel.value(file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'],              checkIfExists: true))
    germline_resource_tbi = Channel.value(file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'],          checkIfExists: true))
    panel_of_normals      = Channel.value(file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'],     checkIfExists: true))
    panel_of_normals_tbi  = Channel.value(file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true))
    interval_file         = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'],              checkIfExists: true))

    BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK ( input, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, interval_file )
}
