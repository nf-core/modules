#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../software/bedtools/complement/main.nf' addParams( options: [:] )
include { BEDTOOLS_GENOMECOV } from '../../../software/bedtools/genomecov/main.nf' addParams( options: [:] )
include { BEDTOOLS_INTERSECT } from '../../../software/bedtools/intersect/main.nf' addParams( options: [:] )
include { BEDTOOLS_MERGE } from '../../../software/bedtools/merge/main.nf' addParams( options: [:] )
include { BEDTOOLS_SLOP} from '../../../software/bedtools/slop/main.nf' addParams( options: [args: '-l 15 -r 30'] )
include { BEDTOOLS_SORT } from '../../../software/bedtools/sort/main.nf' addParams( options: [:] )


workflow test_bedtools_complement {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) ] //metamap

    BEDTOOLS_COMPLEMENT( input )
}

workflow test_bedtools_genomecov {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true) ]

    BEDTOOLS_GENOMECOV( input )
}

workflow test_bedtools_intersect {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/B.bed", checkIfExists: true) ] //metamap

    BEDTOOLS_INTERSECT( input )
}


workflow test_bedtools_merge {
    def input = []
    input = [ [ id:'test' ], // meta map
                file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]
    BEDTOOLS_MERGE(input)
}

workflow test_bedtools_slop {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) ] //metamap
    BEDTOOLS_SLOP ( input )
}

workflow test_bedtools_sort {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]

    BEDTOOLS_SORT( input )

}


