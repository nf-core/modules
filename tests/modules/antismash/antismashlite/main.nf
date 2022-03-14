#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../modules/gunzip/main.nf'
include { ANTISMASHLITEDOWNLOADDATABASES } from '../../../modules/antismash/antismashlitedownloaddatabases/main.nf'
include { ANTISMASHLITE } from '../../../modules/antismash/antismashlite/main.nf'
include { PROKKA } from '../../../modules/prokka/main.nf'

workflow test_antismashlite {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true) ]

    ANTISMASHLITEDOWNLOADDATABASES ()
    ANTISMASHLITE ( input, ANTISMASHLITEDOWNLOADDATABASES.out.db )
}

workflow test_prokka_antismashlite {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]

    GUNZIP ( input )
    PROKKA ( GUNZIP.out.gunzip, [], [] )
    ANTISMASHLITEDOWNLOADDATABASES ()
    ANTISMASH ( PROKKA.out.gbk )
}
