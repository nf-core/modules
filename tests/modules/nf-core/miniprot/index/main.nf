#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MINIPROT_INDEX } from '../../../../../modules/nf-core/miniprot/index/main.nf'

workflow test_miniprot_index {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    MINIPROT_INDEX ( [ [id:'test'], fasta ] )
}
