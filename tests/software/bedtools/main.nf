#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_COMPLEMENT } from '../../../software/bedtools/complement/main.nf' addParams( options: [:] )
include { BEDTOOLS_GENOMECOV } from '../../../software/bedtools/genomecov/main.nf' addParams( options: [:] )
include { BEDTOOLS_INTERSECT } from '../../../software/bedtools/intersect/main.nf' addParams( options: [:] )
include { BEDTOOLS_MERGE } from '../../../software/bedtools/merge/main.nf' addParams( options: [:] )
include { BEDTOOLS_SLOP as BEDTOOLS_SLOP_S} from '../../../software/bedtools/slop/main.nf' addParams( options: [:] )
include { BEDTOOLS_SLOP as BEDTOOLS_SLOP_AS} from '../../../software/bedtools/slop/main.nf' addParams( options: [:] )
include { BEDTOOLS_SORT } from '../../../software/bedtools/sort/main.nf' addParams( options: [publish_dir: 'test_bedtools_sort'] ) // needed for merge


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
              file("${launchDir}/tests/data/bam/test.paired_end.sorted.bam", checkIfExists: true),
              file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) ] //metamap

    BEDTOOLS_GENOMECOV( input )
}

workflow test_bedtools_intersect {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/B.bed", checkIfExists: true) ] //metamap

    BEDTOOLS_INTERSECT( input )
}


//  TODO use output of sort module
workflow test_bedtools_merge {
    def input = []
    input = [ [ id:'test' ], // meta map
                file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]
    BEDTOOLS_MERGE(input)
}

// TODO streamline slop module

// To run with header and pct enabled, type --pct true and --header true with nextflow run command.
/*
Test with l/r method
*/
workflow test_bedtools_slop_asymmetrical {
    def input = []
    input = [ [ id:'test', symmetry: false],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) ] //metamap
    BEDTOOLS_SLOP_AS( input )
}
/*
Test with b method
*/
workflow test_bedtools_slop_symmetrical {
    def input = []
    input = [ [ id:'test', symmetry: true],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true),
              file("${launchDir}/tests/data/bed/genome.sizes", checkIfExists: true) ] //metamap
    BEDTOOLS_SLOP_S( input )
}

workflow test_bedtools_sort {
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bed/A.bed", checkIfExists: true) ]

    BEDTOOLS_SORT( input )

    emit:
    sort = BEDTOOLS_SORT.out.sort

}


