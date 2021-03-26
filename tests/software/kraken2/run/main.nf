#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKEN2_RUN } from '../../../../software/kraken2/run/main.nf' addParams( options: [:] )

workflow test_kraken2_run_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2'], checkIfExists: true)

    KRAKEN2_RUN ( input, db )
}

workflow test_kraken2_run_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    db    = file(params.test_data['sarscov2']['genome']['kraken2'], checkIfExists: true)

    KRAKEN2_RUN ( input, db )
}
