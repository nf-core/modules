#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DESEQ2_DIFFERENTIAL } from '../../../../../modules/nf-core/deseq2/differential/main.nf'

workflow test_deseq2_differential {

    expression_sample_sheet = file(params.test_data['mus_musculus']['genome']['rnaseq_samplesheet'], checkIfExists: true) 
    expression_matrix_file = file(params.test_data['mus_musculus']['genome']['rnaseq_matrix'], checkIfExists: true) 
    expression_contrasts = file(params.test_data['mus_musculus']['genome']['rnaseq_contrasts'], checkIfExists: true) 

    Channel.fromPath(expression_contrasts)
        .splitCsv ( header:true, sep:',' )  
        .map{
            tuple(it, expression_sample_sheet, expression_matrix_file)
        }
        .set{
            input
        }      

    DESEQ2_DIFFERENTIAL (
        input        
    )
}
