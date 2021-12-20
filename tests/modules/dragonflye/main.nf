#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DRAGONFLYE  } from '../../../modules/dragonflye/main.nf'
include { DRAGONFLYE as DRAGONFLYE_RAVEN   } from '../../../modules/dragonflye/main.nf'

workflow test_dragonflye {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/nanopore/subset15000.fq.gz", checkIfExists: true) ]
            ]

    DRAGONFLYE ( input )
}

workflow test_dragonflye_raven {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/nanopore/subset15000.fq.gz", checkIfExists: true) ]
            ]

    DRAGONFLYE_RAVEN ( input )
}
