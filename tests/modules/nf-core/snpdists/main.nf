#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { SNPDISTS } from "$moduleDir/modules/nf-core/snpdists/main.nf"

workflow test_snpdists {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['informative_sites_fas'], checkIfExists: true) ]

    SNPDISTS ( input )
}
