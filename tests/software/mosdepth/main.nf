#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOSDEPTH } from '../../../software/mosdepth/main.nf' addParams( options: [:] )

workflow test_mosdepth {

    input = [ [ id:'test', single_end:true ],
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ],
              [ file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam.bai", checkIfExists: true) ] ]
    dummy = [ file("${launchDir}/tests/data/dummy/dummy_file.txt", checkIfExists: true) ]

    window_size = 100

    MOSDEPTH ( input, dummy, window_size )
}
