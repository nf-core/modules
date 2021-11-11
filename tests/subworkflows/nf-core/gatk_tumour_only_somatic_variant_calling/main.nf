#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_TUMOUR_ONLY_SOMATIC_VARIANT_CALLING } from '../../../../subworkflows/nf-core/gatk_tumour_only_somatic_variant_calling/main' addParams( [:]  )

workflow test_gatk_tumour_only_somatic_variant_calling {
    ch_mutect2_in = [
                    [[ id:'test' ], // meta map
                    [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'],     checkIfExists: true)],
                    [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
                    [] ]
                    ]
    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_fasta'],                            checkIfExists: true)
    fastaidx              = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'],                        checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_dict'],                             checkIfExists: true)
    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz'],                    checkIfExists: true)
    germline_resource_idx = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_vcf_gz_tbi'],                checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz'],           checkIfExists: true)
    panel_of_normals_idx  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_vcf_gz_tbi'],       checkIfExists: true)
    interval_file         = file(params.test_data['homo_sapiens']['genome']['genome_interval_list'],                    checkIfExists: true)

    GATK_TUMOUR_ONLY_SOMATIC_VARIANT_CALLING ( ch_mutect2_in, fasta, fastaidx, dict, germline_resource, germline_resource_idx, panel_of_normals, panel_of_normals_idx, interval_file )
}
