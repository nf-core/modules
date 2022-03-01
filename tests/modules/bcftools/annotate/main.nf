#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_ANNOTATE } from '../../../../modules/bcftools/annotate/main.nf'

workflow test_bcftools_annotate_out_vcf {

    input = [ 
    [ id:'test_compressed_vcf', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    BCFTOOLS_ANNOTATE ( input )
}

workflow test_bcftools_annotate_out_bcf {

    input = [ 
    [ id:'test_compressed_bcf', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_bcf'], checkIfExists: true) ]

    BCFTOOLS_ANNOTATE ( input )
}
