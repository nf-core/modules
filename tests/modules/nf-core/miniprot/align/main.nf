#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { MINIPROT_INDEX } from '../../../../modules/miniprot/index/main.nf'
include { MINIPROT_ALIGN } from '../../../../modules/miniprot/align/main.nf'

workflow test_miniprot_align_gff {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    input_pep = file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)

    MINIPROT_INDEX ( [ [id:'test'], input] )
 
    MINIPROT_ALIGN ( [ [id:'test'], input_pep] , MINIPROT_INDEX.out.index.map { it[1] } )
}

workflow test_miniprot_align_paf {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    input_pep = file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true)

    MINIPROT_INDEX ( [ [id:'test'], input] )
 
    MINIPROT_ALIGN ( [ [id:'test'], input_pep] , MINIPROT_INDEX.out.index.map { it[1] } )
}

