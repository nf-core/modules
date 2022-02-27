#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVDB_MERGE } from '../../../../modules/svdb/merge/main.nf'

workflow test_svdb_merge {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ]
            ]
    priority = [ 'tiddit', 'cnvnator']

    SVDB_MERGE ( input, priority )
}

workflow test_svdb_merge_noprio {

    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test2_vcf'], checkIfExists: true) ]
            ]

    SVDB_MERGE ( input, [] )
}
