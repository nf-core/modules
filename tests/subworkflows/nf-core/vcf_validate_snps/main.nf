#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_VALIDATE_SNPS } from '../../../../subworkflows/nf-core/vcf_validate_snps/main.nf'

workflow test_vcf_validate_snps {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    VCF_VALIDATE_SNPS ( input )
}
