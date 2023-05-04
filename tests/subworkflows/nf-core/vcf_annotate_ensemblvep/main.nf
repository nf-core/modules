#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_DEFAULT } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_CUSTOM  } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'

workflow vcf_annotate_ensemblvep {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)
    ],[
        [ id:'test2' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test2_vcf_gz_tbi'], checkIfExists: true)
    ])

    VCF_ANNOTATE_ENSEMBLVEP_DEFAULT ( input, [[],[]], [[],[]], "WBcel235", "caenorhabditis_elegans", "108", [], [], 5 )
}
