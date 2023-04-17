#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WINDOWMASKER_MKCOUNTS                                } from '../../../../../modules/nf-core/windowmasker/mk_counts/main.nf'
include { WINDOWMASKER_CONVERT as WINDOWMASKER_CONVERT_OASCII  } from '../../../../../modules/nf-core/windowmasker/convert/main.nf'
include { WINDOWMASKER_CONVERT as WINDOWMASKER_CONVERT_BINARY  } from '../../../../../modules/nf-core/windowmasker/convert/main.nf'
include { WINDOWMASKER_CONVERT as WINDOWMASKER_CONVERT_OBINARY } from '../../../../../modules/nf-core/windowmasker/convert/main.nf'

workflow test_windowmasker_convert_oascii {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    WINDOWMASKER_MKCOUNTS ( [ [id:'test'], input ] )
    
    WINDOWMASKER_CONVERT_OASCII ( WINDOWMASKER_MKCOUNTS.out.counts )
}

workflow test_windowmasker_convert_binary {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    WINDOWMASKER_MKCOUNTS ( [ [id:'test'], input ] )
    
    WINDOWMASKER_CONVERT_BINARY ( WINDOWMASKER_MKCOUNTS.out.counts )
}

workflow test_windowmasker_convert_obinary {
    
    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    WINDOWMASKER_MKCOUNTS ( [ [id:'test'], input ] )
    
    WINDOWMASKER_CONVERT_OBINARY ( WINDOWMASKER_MKCOUNTS.out.counts )
}