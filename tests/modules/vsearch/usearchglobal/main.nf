#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_USEARCHGLOBAL } from '../../../../modules/vsearch/usearchglobal/main.nf'

workflow test_vsearch_usearchglobal {
    
    query = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
    db = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    VSEARCH_USEARCHGLOBAL ( query, db, "blast6out_results" )
}
