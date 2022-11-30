#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKV_DOWNLOADDATABASE; CHECKV_DOWNLOADDATABASE as CHECKV_DOWNLOADDATABASE2 } from '../../../../../modules/nf-core/checkv/downloaddatabase/main.nf'

workflow test_checkv_downloaddatabase {
    CHECKV_DOWNLOADDATABASE ([],[])
}

workflow test_checkv_keepdatabase {
    CHECKV_DOWNLOADDATABASE([],[])

    CHECKV_DOWNLOADDATABASE2 ([], CHECKV_DOWNLOADDATABASE.out.checkv_db)

}

workflow test_checkv_updatedatabase {

    input = file(params.test_data['sarscov2']['illumina']['scaffolds_fasta'], checkIfExists: true)

    CHECKV_DOWNLOADDATABASE([],[])

    CHECKV_DOWNLOADDATABASE2( input, CHECKV_DOWNLOADDATABASE.out.checkv_db )

}
