#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UTILS_NFCORE_PIPELINE } from '../../../../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'

workflow test_utils_nfcore_pipeline {
        Channel.of("dummy").collectFile(storeDir: "${params.outdir}/") { it ->
            ["dummy.txt", "it"]
        }
}
