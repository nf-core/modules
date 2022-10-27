#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { RAPIDNJ } from "$moduleDir/modules/nf-core/rapidnj/main.nf"

workflow test_rapidnj {

    input = [ file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    RAPIDNJ ( input )
}
