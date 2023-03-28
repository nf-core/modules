#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INTERPROSCAN } from '../../../../modules/nf-core/interproscan/main.nf'

workflow test_interproscan {
    
    input = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)
    ]
    out_ext = 'tsv'

    INTERPROSCAN ( input, out_ext )
}
