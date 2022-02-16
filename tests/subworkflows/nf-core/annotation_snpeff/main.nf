#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANNOTATION_SNPEFF } from '../../../../subworkflows/nf-core/annotation_snpeff/main'

workflow annotation_snpeff {
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ]

    ANNOTATION_SNPEFF ( input, "WBcel235.99", [] )
}
