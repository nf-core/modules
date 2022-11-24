#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_SNPEFF } from '../../../../subworkflows/nf-core/vcf_annotate_snpeff/main.nf'

workflow vcf_annotate_snpeff {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    VCF_ANNOTATE_SNPEFF ( input, "WBcel235.99", [] )
}
