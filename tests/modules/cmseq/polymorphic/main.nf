#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CMSEQ_POLYMORPHIC } from '../../../../modules/cmseq/polymorphic/main.nf' addParams( options: [:] )

workflow test_cmseq_polymorphic {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    CMSEQ_POLYMORPHIC ( input )
}
