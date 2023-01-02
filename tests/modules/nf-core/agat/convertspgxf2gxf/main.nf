#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_CONVERTSPGXF2GXF } from '../../../../../modules/nf-core/agat/convertspgxf2gxf/main.nf'

workflow test_gff {

    input = [
        [ id: 'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) 
    ]

    AGAT_CONVERTSPGXF2GXF ( input )
}

workflow test_gtf {

    input = [
        [ id: 'test' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) 
    ]

    AGAT_CONVERTSPGXF2GXF ( input )
}
