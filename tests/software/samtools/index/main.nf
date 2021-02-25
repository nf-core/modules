#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../software/samtools/index/main.nf' addParams( options: [:] )

workflow test_samtools_index {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
    SAMTOOLS_INDEX ( input )
}
