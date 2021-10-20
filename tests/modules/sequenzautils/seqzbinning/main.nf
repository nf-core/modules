#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQUENZAUTILS_SEQZBINNING } from '../../../../modules/sequenzautils/seqzbinning/main.nf' addParams( options: [:] )

workflow test_sequenzautils_seqzbinning {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_seqz'], checkIfExists: true) ]

    SEQUENZAUTILS_SEQZBINNING ( input )
}
