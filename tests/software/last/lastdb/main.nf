#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LAST_LASTDB } from '../../../../software/last/lastdb/main.nf' addParams( options: ['args': '-Q0'] )

workflow test_last_lastdb {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
            ]

    LAST_LASTDB ( input )
}

workflow test_last_lastdb_gzipped_input {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]

    LAST_LASTDB ( input )
}
