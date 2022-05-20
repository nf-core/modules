#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_RHOCALL } from '../../../../modules/bcftools/rhocall/main.nf'

workflow test_bcftools_rhocall {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BCFTOOLS_RHOCALL ( input )
}
