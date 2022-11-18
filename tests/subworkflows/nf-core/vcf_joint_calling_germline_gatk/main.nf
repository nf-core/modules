#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_JOINT_CALLING_GERMLINE_GATK } from '../../../../subworkflows/nf-core/vcf_joint_calling_germline_gatk/main.nf'

workflow test_vcf_joint_calling_germline_gatk {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    VCF_JOINT_CALLING_GERMLINE_GATK ( input )
}
