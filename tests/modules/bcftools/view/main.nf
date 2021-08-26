#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_VIEW } from '../../../../modules/bcftools/view/main.nf' addParams( options: ['args': '--no-version'] )

workflow test_bcftools_view {

    bed = []
    targets = []
    samples = []

    input = [[ id:'test2', single_end:false ], // meta map
             [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)]]

    BCFTOOLS_VIEW ( input, bed, targets, samples )
}

workflow test_bcftools_view_with_optional_files {

    bed = []
    targets = []
    samples = []

    input = [[ id:'test2', single_end:false ], // meta map
             [ file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true)]]

    BCFTOOLS_VIEW ( input, bed, targets, samples )
}
