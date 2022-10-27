#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { SEGEMEHL_INDEX } from "$moduleDir/modules/nf-core/segemehl/index/main.nf"

workflow test_segemehl_index {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    SEGEMEHL_INDEX ( fasta )
}
