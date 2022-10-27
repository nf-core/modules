#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { BBMAP_INDEX } from "$moduleDir/modules/nf-core/bbmap/index/main.nf"

workflow test_bbmap_index {

    input = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BBMAP_INDEX ( input )
}
