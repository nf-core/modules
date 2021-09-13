#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAPHOPIASCCMEC } from '../../../modules/staphopiasccmec/main.nf' addParams( options: [:] )

workflow test_staphopiasccmec {
    
    input = [ [ id:'test' ],
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    STAPHOPIASCCMEC ( input )
}
