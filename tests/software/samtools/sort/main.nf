#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_SORT } from '../../../../software/samtools/sort/main.nf' addParams( options: [:] )

workflow test_samtools_sort  {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/bam/sarscov2_paired_aln.bam", checkIfExists: true) ]
    SAMTOOLS_SORT ( input )
}
