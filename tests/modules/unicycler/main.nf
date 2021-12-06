#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNICYCLER } from '../../../modules/unicycler/main.nf'

workflow test_unicycler_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true) ],
              []
            ]

    UNICYCLER ( input )
}

workflow test_unicycler_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['bacteroides_fragilis']['illumina']['test1_2_fastq_gz'], checkIfExists: true) ],
              []
            ]

    UNICYCLER ( input )
}

workflow test_unicycler_shortreads_longreads {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['bacteroides_fragilis']['illumina']['test1_2_fastq_gz'], checkIfExists: true) ],
              [ file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    UNICYCLER ( input )
}
