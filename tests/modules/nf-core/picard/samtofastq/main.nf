#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SAMTOFASTQ } from '../../../../../modules/nf-core/picard/samtofastq/main.nf'

workflow test_picard_samtofastq_single {
    input = [ [ id:'test', single_end: true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
            ]
    PICARD_SAMTOFASTQ ( input )
}

workflow test_picard_samtofastq_paired {
    input = [ [ id:'test', single_end: false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
            ]
    PICARD_SAMTOFASTQ ( input )
}

workflow test_picard_samtofastq_paired_keep_unpaired {
    input = [ [ id:'test', single_end: false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
            ]
    PICARD_SAMTOFASTQ ( input )
}
