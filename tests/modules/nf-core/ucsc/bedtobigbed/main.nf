#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { UCSC_BEDTOBIGBED } from "$moduleDir/modules/nf-core/ucsc/bedtobigbed/main.nf"
include { UCSC_BEDTOBIGBED as UCSC_BEDTOBIGBED_AS } from "$moduleDir/modules/nf-core/ucsc/bedtobigbed/main.nf"

workflow test_ucsc_bedtobigbed {
    input = [ [ id: 'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true ) ]
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_BEDTOBIGBED ( input, sizes, [] )
}

workflow test_ucsc_bedtobigbed_autosql {
    input = [ [ id: 'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true ) ]
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    autosql = file(params.test_data['sarscov2']['genome']['bed_autosql'], checkIfExists: true)

    UCSC_BEDTOBIGBED_AS ( input, sizes, autosql )
}
