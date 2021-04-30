#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNICYCLER } from '../../../software/unicycler/main.nf' addParams( options: [:] )

workflow test_unicycler_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    UNICYCLER ( input )
}

workflow test_unicycler_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    UNICYCLER ( input )
}
