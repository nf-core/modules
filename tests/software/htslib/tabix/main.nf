#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HTSLIB_TABIX as TABIX_BED } from '../../../../software/htslib/tabix/main.nf' addParams( options: ['args': '-p bed'] )
include { HTSLIB_TABIX as TABIX_GFF } from '../../../../software/htslib/tabix/main.nf' addParams( options: ['args': '-p gff'] )
include { HTSLIB_TABIX as TABIX_VCF } from '../../../../software/htslib/tabix/main.nf' addParams( options: ['args': '-p vcf'] )

workflow test_htslib_tabix_bed {
    TABIX_BED ( file("${launchDir}/tests/data/bed/B.bed.gz", checkIfExists: true) )
}

workflow test_htslib_tabix_gff {
    TABIX_GFF ( file("${launchDir}/tests/data/gff/a.gff3.gz", checkIfExists: true) )
}

workflow test_htslib_tabix_vcf {
    TABIX_VCF ( file("${launchDir}/tests/data/vcf/test.vcf.gz", checkIfExists: true) )
}
