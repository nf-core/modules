#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV_UPDATEDATABASE }   from '../../../../../modules/nf-core/checkv/updatedatabase/main.nf'
include { CHECKV_DOWNLOADDATABASE } from '../../../../../modules/nf-core/checkv/downloaddatabase/main.nf'

workflow test_checkv_updatedatabase {



    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true)]

    CHECKV_DOWNLOADDATABASE( input[0])

    CHECKV_UPDATEDATABASE ( input ,CHECKV_DOWNLOADDATABASE.out.checkv_db )
}
