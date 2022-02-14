#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_ANNOTATE } from '../../../../modules/bcftools/annotate/main.nf'

workflow test_bcftools_annotate {
    
    input = [ 
    [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true) ]

    BCFTOOLS_ANNOTATE ( input )
}
