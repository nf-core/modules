#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UTILS_NEXTFLOW_PIPELINE } from '../../../../subworkflows/nf-core/utils_nextflow_pipeline/main.nf'

workflow test_utils_nextflow_pipeline {
        Channel.of("dummy").collectFile(storeDir: "${params.outdir}/") { it ->
            ["dummy.txt", "it"]
        }
}
