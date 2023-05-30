#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IPHOP_ADDTODB } from '../../../../../modules/nf-core/iphop/addtodb/main.nf'

workflow test_iphop_addtodb {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    IPHOP_ADDTODB ( input )
}
