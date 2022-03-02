#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MTNUCRATIO } from '../../../modules/mtnucratio/main.nf'

workflow test_mtnucratio {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true)]
    mt_id = 'mt_id'

    MTNUCRATIO ( input, mt_id )
}
