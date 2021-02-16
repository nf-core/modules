#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HTSLIB_TABIX } from '../../../../software/htslib/tabix/main.nf' addParams( options: [:] )

workflow test_htslib_tabix {
    HTSLIB_TABIX ( file("${launchDir}/tests/data/vcf/test.vcf.gz", checkIfExists: true) )
}