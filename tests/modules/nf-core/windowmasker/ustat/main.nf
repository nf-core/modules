#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { WINDOWMASKER_MKCOUNTS } from '../../../../../modules/nf-core/windowmasker/mk_counts/main.nf'
include { WINDOWMASKER_USTAT    } from '../../../../../modules/nf-core/windowmasker/ustat/main.nf'


workflow test_windowmasker_ustat {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    WINDOWMASKER_MKCOUNTS ( [ [id:'test'], input ] )

    WINDOWMASKER_USTAT ( WINDOWMASKER_MKCOUNTS.out.counts, [ [id:'test'], input ])

}