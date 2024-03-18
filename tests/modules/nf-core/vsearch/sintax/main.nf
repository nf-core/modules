#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_SINTAX } from '../../../../../modules/nf-core/vsearch/sintax/main.nf'

workflow test_vsearch_sintax {
    
    query = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    db = file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)

    VSEARCH_SINTAX ( [[id:'test'], query], db )
}
