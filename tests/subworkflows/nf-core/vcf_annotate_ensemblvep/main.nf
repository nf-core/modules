#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_DEFAULT } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_CUSTOM  } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'

workflow vcf_annotate_ensemblvep {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true), []
    ])

    VCF_ANNOTATE_ENSEMBLVEP_DEFAULT ( input, [[],[]], "WBcel235", "caenorhabditis_elegans", "108", [], [] )
}

workflow vcf_annotate_ensemblvep_custom {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
        [
            file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test3_vcf'], checkIfExists: true)
        ]
    ])

    VCF_ANNOTATE_ENSEMBLVEP_CUSTOM ( input, [[],[]], "WBcel235", "caenorhabditis_elegans", "108", [], [] )
}
