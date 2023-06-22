#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMOGRAPH } from '../../../../modules/nf-core/chromograph/main.nf'

workflow test_chromograph {

    sites = [
        [ id:'test', single_end:false ], // meta map
        file("/home/ramprasad.neethiraj/nextflow/test-datasets/data/genomics/homo_sapiens/genome/updsites.bed", checkIfExists: true)
    ]

    CHROMOGRAPH ( [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], sites )
}
