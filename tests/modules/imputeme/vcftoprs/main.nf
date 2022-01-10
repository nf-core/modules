#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IMPUTEME_VCFTOPRS } from '../../../../modules/imputeme/vcftoprs/main.nf'

workflow test_imputeme_vcftoprs {
    
    input = [ 
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true)
    ]

    IMPUTEME_VCFTOPRS ( input )
}
