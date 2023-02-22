#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEROBA_RUN } from '../../../../../modules/nf-core/seroba/run/main.nf'

workflow test_seroba_run {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/ERR044595_1M_1.fastq.gz", checkIfExists: true),
                file("https://github.com/nf-core/test-datasets/raw/bacass/ERR044595_1M_2.fastq.gz", checkIfExists: true) ]
            ]

    SEROBA_RUN ( input )
}
