#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANTISMASH } from '../../../modules/antismash/main.nf' addParams( options: [:] )

workflow test_antismash {

    input = [ [ id:'test' ], // meta map
            file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]
    ANTISMASH ( input )
}
