#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEQC_READDUPLICATION }   from '../../../../modules/rseqc/readduplication/main.nf'

workflow test_rseqc_readduplication {
    input = [ [ id:'test', single_end: false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
            ]

    RSEQC_READDUPLICATION ( input )
}
