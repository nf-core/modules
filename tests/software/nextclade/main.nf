#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NEXTCLADE } from '../../../software/nextclade/main.nf' addParams( options: [:] )

workflow test_nextclade_json {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    NEXTCLADE ( input, 'json' )
}

workflow test_nextclade_csv {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    NEXTCLADE ( input, 'csv' )
}

workflow test_nextclade_tsv {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    NEXTCLADE ( input, 'tsv' )
}

workflow test_nextclade_tree {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]
    NEXTCLADE ( input, 'tree' )
}
