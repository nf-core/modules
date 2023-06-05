#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGAT_SPSTATISTICS } from '../../../../../modules/nf-core/agat/spstatistics/main.nf'

workflow test_statistics {

    gff = [ file(params.test_data['sarscov2']['genome']['genome_gff3'], checkIfExists: true) ]


    AGAT_SPSTATISTICS ( [ [id:'test'], gff ] )

}

