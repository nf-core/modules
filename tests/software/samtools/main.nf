#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from '../../../software/samtools/flagstat/main.nf' addParams( options: [:] )
include { SAMTOOLS_IDXSTATS } from '../../../software/samtools/idxstats/main.nf' addParams( options: [:] )
include { SAMTOOLS_INDEX } from '../../../software/samtools/index/main.nf' addParams( options: [:] )
include { SAMTOOLS_SORT } from '../../../software/samtools/sort/main.nf' addParams( options: [:] )
include { SAMTOOLS_STATS } from '../../../software/samtools/stats/main.nf' addParams( options: [:] )
include { SAMTOOLS_VIEW } from '../../../software/samtools/view/main.nf' addParams( options: [:] )
include { SAMTOOLS_MPILEUP } from '../../../software/samtools/mpileup/main.nf' addParams( options: [:] )
include { SAMTOOLS_FAIDX } from '../../../software/samtools/faidx/main.nf' addParams( options: [:] )

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

// FIXME Why is this passing it an already sorted bam?
workflow test_samtools_sort  {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]

    SAMTOOLS_SORT ( input )
}

workflow test_samtools_stats {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam.bai", checkIfExists: true) ]

    SAMTOOLS_STATS ( input )
}

workflow test_samtools_view {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]

    SAMTOOLS_VIEW ( input )
}

workflow test_samtools_mpileup {

    def input = []
    def fasta = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]
    fasta = [  file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)  ]

    SAMTOOLS_MPILEUP ( input, fasta )
}

workflow test_samtools_faidx {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true) ]

    SAMTOOLS_FAIDX ( input )
}
