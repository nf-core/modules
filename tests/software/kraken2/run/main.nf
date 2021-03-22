#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKEN2_RUN } from '../../../../software/kraken2/run/main.nf' addParams( options: [:] )

workflow test_kraken2_run_single_end {
    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
                file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_1.fastq.gz", checkIfExists: true) ]

    db = [ file("${launchDir}/tests/data/genomics/sarscov2/db/kraken2", checkIfExists: true) ]

    KRAKEN2_RUN ( input, db )
}

workflow test_kraken2_run_paired_end {
    def input, db = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/fastq/test_*.fastq.gz", checkIfExists: true) ]

    db = [ file("${launchDir}/tests/data/genomics/sarscov2/db/kraken2", checkIfExists: true) ]

    KRAKEN2_RUN ( input, db )
}
