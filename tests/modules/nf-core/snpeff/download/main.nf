#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPEFF_DOWNLOAD } from '../../../../../modules/nf-core/snpeff/download/main.nf'

workflow test_snpeff_download {
    input = [ [ id:'test' ], "WBcel235", "105"]

    SNPEFF_DOWNLOAD ( input )
}
