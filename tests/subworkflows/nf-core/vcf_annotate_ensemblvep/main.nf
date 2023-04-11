#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'

workflow vcf_annotation_ensemblvep {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ])

    custom_extra_files = Channel.of([
        [ id:'test' ],
        [
            file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test3_vcf'], checkIfExists: true)            
        ]
    ])

    VCF_ANNOTATE_ENSEMBLVEP ( input, [[],[]], "WBcel235", "caenorhabditis_elegans", "108", [[],[]], [], custom_extra_files )
}
