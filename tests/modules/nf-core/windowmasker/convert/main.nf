#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WINDOWMASKER_MKCOUNTS } from '../../../../../modules/nf-core/windowmasker/mk_counts/main.nf'
include { WINDOWMASKER_CONVERT  } from '../../../../../modules/nf-core/windowmasker/convert/main.nf'

workflow test_windowmasker_convert {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    WINDOWMASKER_MKCOUNTS ( [ [id:'test'], input ] )
    
    WINDOWMASKER_CONVERT ( WINDOWMASKER_MKCOUNTS.out.counts )
}
