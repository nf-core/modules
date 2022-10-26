#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_STATISTICS } from '../../../../../modules/nf-core/agat/statistics/main.nf'

workflow test_statistics {

    gff = [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) ]


    AGAT_STATISTICS ( [ [id:'test'], gff ] )

}

workflow test_basic {

    gff = [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) ]


    AGAT_STATISTICS ( [ [id:'test'], gff ] )

}

