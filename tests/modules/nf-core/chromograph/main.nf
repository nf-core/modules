#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHROMOGRAPH as CHROMOGRAPH_SITES   } from '../../../../modules/nf-core/chromograph/main.nf'
include { CHROMOGRAPH as CHROMOGRAPH_REGIONS } from '../../../../modules/nf-core/chromograph/main.nf'

workflow test_chromograph_sites {

    sites = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['updsites_bed'], checkIfExists: true)
    ]

    CHROMOGRAPH_SITES ( [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], sites )
}

workflow test_chromograph_regions {

    regions = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['updsites_bed'], checkIfExists: true)
    ]

    CHROMOGRAPH_REGIONS ( [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], regions, [[],[]] )
}
