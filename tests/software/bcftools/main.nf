#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//keep arg of bcftools_filter, otherwise md5 will change on each execution
include { BCFTOOLS_FILTER } from '../../../software/bcftools/filter/main.nf' addParams( options: ['args': '--no-version'] )
include { BCFTOOLS_STATS } from '../../../software/bcftools/stats/main.nf' addParams( options: [:] )
include { BCFTOOLS_BGZIP } from '../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX } from '../../../software/bcftools/tabix/main.nf' addParams( options: [:] )
include { BCFTOOLS_CONSENSUS } from '../../../software/bcftools/consensus/main.nf' addParams( options: [:] )
include { BCFTOOLS_ISEC } from '../../../software/bcftools/isec/main.nf' addParams( options: ['args': '--nfiles +2 --output-type z --no-version'] )
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP2 } from '../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX as BCFTOOLS_TABIX2 } from '../../../software/bcftools/tabix/main.nf' addParams( options: [:] )
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP3 } from '../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX as BCFTOOLS_TABIX3 } from '../../../software/bcftools/tabix/main.nf' addParams( options: [:] )

workflow test_bcftools_filter {

    def input = []
    input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    BCFTOOLS_FILTER ( input )
}

workflow test_bcftools_stats {

    def input = []
    input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    BCFTOOLS_STATS ( input )
}

workflow test_bcftools_bgzip {

    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    BCFTOOLS_BGZIP ( input )
}

workflow test_bcftools_tabix {

    def input = []
    input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    BCFTOOLS_BGZIP ( input )
    BCFTOOLS_TABIX ( BCFTOOLS_BGZIP.out.vcf )
}

workflow test_bcftools_consensus {

    def input = []
    input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]
    def fasta = []
    fasta = Channel.of ([ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/vcf/test.consensus.fa", checkIfExists: true) ])

    BCFTOOLS_BGZIP ( input )
    BCFTOOLS_TABIX ( BCFTOOLS_BGZIP.out.vcf )
    BCFTOOLS_CONSENSUS ( BCFTOOLS_BGZIP.out.vcf
                                    .join( BCFTOOLS_TABIX.out.tbi, by: [0] )
                                    .join( fasta, by: [0] ) )
}

workflow test_bcftools_isec {

    def input = []
    input1 = [ [ id:'test1' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    input2 = [ [ id:'test2' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    input3 = [ [ id:'test3' ], // meta map
            file("${launchDir}/tests/data/vcf/test.vcf", checkIfExists: true) ]

    BCFTOOLS_BGZIP  ( input1 )
    BCFTOOLS_TABIX  ( BCFTOOLS_BGZIP.out.vcf )
    BCFTOOLS_BGZIP2 ( input2 )
    BCFTOOLS_TABIX2 ( BCFTOOLS_BGZIP2.out.vcf )
    BCFTOOLS_BGZIP3 ( input3 )
    BCFTOOLS_TABIX3 ( BCFTOOLS_BGZIP3.out.vcf )
    BCFTOOLS_ISEC   ( BCFTOOLS_BGZIP.out.vcf
                        .mix( BCFTOOLS_TABIX.out.tbi )
                        .mix( BCFTOOLS_BGZIP2.out.vcf )
                        .mix( BCFTOOLS_TABIX2.out.tbi )
                        .mix( BCFTOOLS_BGZIP3.out.vcf )
                        .mix( BCFTOOLS_TABIX3.out.tbi)
                        .map { [ it[1] ]}
                        .collect()
                        .map { [ [id: 'test'], it ] } )
}
