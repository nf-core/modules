#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_CONVERTSPGFF2GTF } from '../../../../../modules/nf-core/agat/convertspgff2gtf/main.nf'

workflow test_gff {

    input = [ 
        [ id: 'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) 
    ]

    AGAT_CONVERTSPGFF2GTF ( input )
}
