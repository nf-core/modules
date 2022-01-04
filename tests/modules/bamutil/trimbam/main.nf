#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMUTIL_TRIMBAM } from '../../../../modules/bamutil/trimbam/main.nf'

workflow test_bamutil_trimbam {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
              2,
              2 ]

    BAMUTIL_TRIMBAM ( input )
}
