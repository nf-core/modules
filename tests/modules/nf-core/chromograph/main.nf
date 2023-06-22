#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMOGRAPH } from '../../../../modules/nf-core/chromograph/main.nf'

workflow test_chromograph {

    sites = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['updsites_bed'], checkIfExists: true)
    ]

    CHROMOGRAPH ( [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], sites )
}
