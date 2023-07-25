#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPEFF_DOWNLOAD } from '../../../../../modules/nf-core/snpeff/download/main.nf'

workflow test_snpeff_download {
    snpeff_cache_version = "110"
    snpeff_genome = "WBcel235"
    snpeff_cache_input = Channel.of([[id:"${snpeff_genome}.${snpeff_cache_version}"], snpeff_genome, snpeff_cache_version])

    SNPEFF_DOWNLOAD(snpeff_cache_input)
}
