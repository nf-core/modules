#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TABIX_TABIX as TABIX_BED } from '../../../../software/tabix/tabix/main.nf' addParams( options: ['args': '-p bed'] )
include { TABIX_TABIX as TABIX_GFF } from '../../../../software/tabix/tabix/main.nf' addParams( options: ['args': '-p gff'] )
include { TABIX_TABIX as TABIX_VCF } from '../../../../software/tabix/tabix/main.nf' addParams( options: ['args': '-p vcf'] )

workflow test_tabix_tabix_bed {
    def input = []
    input = [ [ id:'B.bed' ], // meta map
              [ file("${launchDir}/tests/data/genomics/sarscov2/bed/sarscov2.bed.gz", checkIfExists: true) ] ]

    TABIX_BED ( input )
}

workflow test_tabix_tabix_gff {
    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/generic/a.gff3.gz", checkIfExists: true) ] ]

    TABIX_GFF ( input )
}

workflow test_tabix_tabix_vcf {
    def input = []
    input = [ [ id:'test.vcf' ], // meta map
              [ file("${launchDir}/tests/data/generic/vcf/test.vcf.gz", checkIfExists: true) ] ]

    TABIX_VCF ( input )
}
