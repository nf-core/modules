#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV_DOWNLOADDATABASE } from '../../../../../modules/nf-core/checkv/downloaddatabase/main.nf'

workflow test_checkv_downloaddatabase {
    CHECKV_DOWNLOADDATABASE ([],[])

    emit:
    checkv_db = CHECKV_DOWNLOADDATABASE.out.checkv_db // path: checkv_db/
}

workflow test_checkv_keepdatabase {

    CHECKV_DOWNLOADDATABASE ([], test_checkv_downloaddatabase.out.checkv_db)

}

workflow test_checkv_updatedatabase {

    input = file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)

    CHECKV_DOWNLOADDATABASE ( input, test_checkv_downloaddatabase.out.checkv_db )

}
