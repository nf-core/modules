#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV_ENDTOEND }         from '../../../../../modules/nf-core/checkv/endtoend/main.nf'
include { CHECKV_DOWNLOADDATABASE } from '../../../../../modules/nf-core/checkv/downloaddatabase/main.nf'


workflow test_checkv_endtoend {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true) ]

    CHECKV_DOWNLOADDATABASE ()

    CHECKV_ENDTOEND( input, CHECKV_DOWNLOADDATABASE.out.checkv_db )
}
