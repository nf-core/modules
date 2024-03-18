#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_MAKEBINS } from '../../../../../modules/nf-core/cooler/makebins/main.nf'

workflow test_cooler_makebins {

    input = [ [ id:'test' ], // meta map
        file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true),
        "1000000"
    ]

    COOLER_MAKEBINS ( input )
}
