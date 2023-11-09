#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DRAGONFLYE  } from '../../../../modules/nf-core/dragonflye/main.nf'
include { DRAGONFLYE as DRAGONFLYE_RAVEN   } from '../../../../modules/nf-core/dragonflye/main.nf'
include { DRAGONFLYE as DRAGONFLYE_HYBRID   } from '../../../../modules/nf-core/dragonflye/main.nf'

workflow test_dragonflye {
    input = [ [ id:'test', single_end:true ], // meta map
              [], // short reads for polishing
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/nanopore/subset15000.fq.gz", // long reads
              checkIfExists: true) ]
            ]

    DRAGONFLYE ( input )
}

workflow test_dragonflye_raven {
    input = [ [ id:'test', single_end:true ], // meta map
              [], // short reads for polishing
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/nanopore/subset15000.fq.gz", // long reads
              checkIfExists: true) ]
            ]

    DRAGONFLYE_RAVEN ( input )
}

workflow test_dragonflye_hybrid {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/A1403KPN_R1.fastq.gz"), file("https://github.com/nf-core/test-datasets/raw/bacass/A1403KPN_R2.fastq.gz") ], // short reads for polishing
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/nanopore/A1403KPN.fq.gz", // long reads
              checkIfExists: true) ]
            ]

    DRAGONFLYE_HYBRID ( input )
}
