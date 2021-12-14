#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAPHOPIASCCMEC } from '../../../modules/staphopiasccmec/main.nf'
include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_HAMMING } from '../../../modules/staphopiasccmec/main.nf'

workflow test_staphopiasccmec {
    
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    STAPHOPIASCCMEC ( input )
}

workflow test_staphopiasccmec_hamming {
    
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    STAPHOPIASCCMEC_HAMMING ( input )
}
