#!/usr/bin/env nextflow



include { CUSTOM_GETCHROMSIZES } from '../../../../modules/custom/getchromsizes/main.nf'

workflow test_custom_getchromsizes {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 

    CUSTOM_GETCHROMSIZES ( input )
}
