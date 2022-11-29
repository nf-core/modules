#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_SQSTATBASIC } from '../../../../../modules/nf-core/agat/sqstatbasic/main.nf'

workflow test_basic {

    gff = [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) ]


    AGAT_SQSTATBASIC ( [ [id:'test'], gff ] )

}

