#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROTEINORTHO } from '../../../../modules/nf-core/proteinortho/main.nf' addParams( options: [args: '-e=10'] )

workflow test_proteinortho {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ 
            file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true),
            file(params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true)
        ]
    ]

    PROTEINORTHO ( input )
}
