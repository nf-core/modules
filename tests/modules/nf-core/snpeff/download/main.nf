#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPEFF_DOWNLOAD } from '../../../../../modules/nf-core/snpeff/download/main.nf'

snpeff_cache_version = "105"
snpeff_genome = "WBcel235"
snpeff_cache_input = Channel.of([[id:"${snpeff_genome}.${snpeff_cache_version}"], snpeff_genome, snpeff_cache_version])

workflow test_snpeff_download {
    SNPEFF_DOWNLOAD(snpeff_cache_input)
}
