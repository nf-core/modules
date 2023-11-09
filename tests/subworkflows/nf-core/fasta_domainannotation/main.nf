#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../../../../modules/nf-core/diamond/makedb/main.nf'
include { FASTA_DOMAINANNOTATION } from '../../../../subworkflows/nf-core/fasta_domainannotation/main.nf'

workflow test_fasta_domainannotation_diamond {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    input = Channel.of( [ [id:'test'], fasta ] )
    blast_fasta = Channel.of( [ [id:'test'], fasta ] )
    blast_mode = "diamond"

    FASTA_DOMAINANNOTATION ( input, blast_fasta, blast_mode )
}

workflow test_fasta_domainannotation_blast {

    fasta = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]
    input = Channel.of( [ [id:'test'], fasta ] )
    blast_fasta = Channel.of( [ [id:'test'], fasta ] )
    blast_mode = "blast"

    FASTA_DOMAINANNOTATION ( input, blast_fasta, blast_mode )
}
