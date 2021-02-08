#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_BGZIP } from '../../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX } from '../../../../software/bcftools/tabix/main.nf' addParams( options: [:] )
include { BCFTOOLS_ISEC } from '../../../../software/bcftools/isec/main.nf' addParams( options: ['args': '--nfiles +2 --output-type z --no-version'] )
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP2 } from '../../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX as BCFTOOLS_TABIX2 } from '../../../../software/bcftools/tabix/main.nf' addParams( options: [:] )
include { BCFTOOLS_BGZIP as BCFTOOLS_BGZIP3 } from '../../../../software/bcftools/bgzip/main.nf' addParams( options: [:] )
include { BCFTOOLS_TABIX as BCFTOOLS_TABIX3 } from '../../../../software/bcftools/tabix/main.nf' addParams( options: [:] )

workflow test_bcftools_isec {

    def input1, input2, input3 = []

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

    vcfs = BCFTOOLS_BGZIP.out.vcf
                        .mix( BCFTOOLS_BGZIP2.out.vcf )
                        .mix( BCFTOOLS_BGZIP3.out.vcf )
                        .map { [ it[1] ]}.collect().map { [ [id: 'test'], it ] }

    tbis = BCFTOOLS_TABIX.out.tbi
                        .mix( BCFTOOLS_TABIX2.out.tbi )
                        .mix( BCFTOOLS_TABIX3.out.tbi )
                        .map { [ it[1] ]}.collect().map { [ [id: 'test'], it ] }

    BCFTOOLS_ISEC ( vcfs.join(tbis) )
}
