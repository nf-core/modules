#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CMSEQ_POLYMUT } from '../../../../modules/cmseq/polymut/main.nf' addParams( options: [:] )

workflow test_cmseq_polymut{
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    CMSEQ_POLYMUT( input )
}
