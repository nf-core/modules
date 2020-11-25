#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from '../../../software/samtools/flagstat/main.nf' addParams( options: [:] )
include { SAMTOOLS_IDXSTATS } from '../../../software/samtools/idxstats/main.nf' addParams( options: [:] )
include { SAMTOOLS_INDEX } from '../../../software/samtools/index/main.nf' addParams( options: [:] )

workflow test_samtools_flagstat {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam.bai", checkIfExists: true) ]

    SAMTOOLS_FLAGSTAT ( input )
}

workflow test_samtools_idxstats {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam.bai", checkIfExists: true) ]

    SAMTOOLS_IDXSTATS ( input )
}

workflow test_samtools_index {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]

    SAMTOOLS_INDEX ( input )
}
