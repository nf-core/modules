#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_IDXSTATS } from '../../../../software/samtools/idxstats/main.nf' addParams( options: [:] )

workflow test_samtools_idxstats {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test-sc2-artic-v3-sorted-trimmed.bam.bai", checkIfExists: true) ]
    SAMTOOLS_IDXSTATS ( input )
}
