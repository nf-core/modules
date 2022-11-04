#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIPROT_INDEX } from '../../../../modules/miniprot/index/main.nf'

workflow test_miniprot_index {
    fasta =file(params.tol_test_data['small_genome']['Oscheius_sp']['assembly']['assembly_fasta'],checkIfExists: true)
    MINIPROT_INDEX ( [ [id:'test'], fasta ] )
}
