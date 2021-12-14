#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../modules/gunzip/main.nf'
include { ANTISMASH } from '../../../modules/antismash/main.nf'
include { PROKKA } from '../../../modules/prokka/main.nf'

workflow test_antismash {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true) ]

    ANTISMASH ( input )
}

workflow test_prokka_antismash {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]

    GUNZIP ( input )
    PROKKA ( GUNZIP.out.gunzip, [], [] )
    ANTISMASH ( PROKKA.out.gbk )
}