#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MEDAKA } from '../../../modules/medaka/main.nf'

workflow test_medaka {
    
    input = [ 
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    
    MEDAKA ( input )
}
