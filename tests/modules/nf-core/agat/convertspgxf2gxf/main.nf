#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_CONVERTSPGXF2GXF } from '../../../../../modules/nf-core/agat/convertspgxf2gxf/main.nf'

workflow test_gff {

    gff = [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) ]

    AGAT_CONVERTSPGXF2GXF ( [ [id:'test'], gff ] )
}

workflow test_gtf {

    gff = [ file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true) ]

    AGAT_CONVERTSPGXF2GXF ( [ [id:'test'], gff ] )
}
