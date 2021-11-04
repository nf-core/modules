#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_GETCHROMSIZES } from '../../../../modules/samtools/getchromsizes/main.nf' addParams( options: [:] )

workflow test_samtools_getchromsizes {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) 

    SAMTOOLS_GETCHROMSIZES ( input )
}
