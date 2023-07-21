#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_DEFAULT } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_CUSTOM  } from '../../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { ENSEMBLVEP_DOWNLOAD                                        } from '../../../../../modules/nf-core/ensemblvep/download/main.nf'

workflow vcf_annotate_ensemblvep {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true), []
    ])

    input_vep_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_vep_cache)

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    VCF_ANNOTATE_ENSEMBLVEP_DEFAULT ( input, [[],[]], "WBcel235", "caenorhabditis_elegans", "110", cache, [] )
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

    input_cache = [[id:"test"], "WBcel235", "caenorhabditis_elegans", "110"]

    ENSEMBLVEP_DOWNLOAD(input_cache)

    cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }

    VCF_ANNOTATE_ENSEMBLVEP_CUSTOM ( input, [[],[]], "WBcel235", "caenorhabditis_elegans", "110", cache, [] )
}
