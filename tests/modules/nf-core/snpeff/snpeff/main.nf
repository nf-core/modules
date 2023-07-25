#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNPEFF_DOWNLOAD } from '../../../../../modules/nf-core/snpeff/download/main'
include { SNPEFF_SNPEFF   } from '../../../../../modules/nf-core/snpeff/snpeff/main'

snpeff_cache_version = "110"
snpeff_genome = "WBcel235"
snpeff_cache_input = Channel.of([[id:"${snpeff_genome}.${snpeff_cache_version}"], snpeff_genome, snpeff_cache_version])

workflow test_snpeff_snpeff {
    input = Channel.of([
        [ id:'test' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
    ])

    SNPEFF_DOWNLOAD(snpeff_cache_input)

    snpeff_cache = SNPEFF_DOWNLOAD.out.cache.first()

    SNPEFF_SNPEFF(input, "${snpeff_genome}.${snpeff_cache_version}", snpeff_cache)
}
