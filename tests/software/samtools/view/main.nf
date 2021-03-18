#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_VIEW } from '../../../../software/samtools/view/main.nf' addParams( options: [:] )

workflow test_samtools_view {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]
    SAMTOOLS_VIEW ( input )
}
