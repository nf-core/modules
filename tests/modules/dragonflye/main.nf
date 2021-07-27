#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DRAGONFLYE                       } from '../../../modules/dragonflye/main.nf' addParams( options: [args: '--assembler flye --gsize 30000']    )
include { DRAGONFLYE as DRAGONFLYE_MINIASM } from '../../../modules/dragonflye/main.nf' addParams( options: [args: '--assembler miniasm --gsize 30000'] )
include { DRAGONFLYE as DRAGONFLYE_RAVEN   } from '../../../modules/dragonflye/main.nf' addParams( options: [args: '--assembler raven --gsize 30000']   )

workflow test_dragonflye {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    DRAGONFLYE ( input )
}

workflow test_dragonflye_megahit {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    DRAGONFLYE_MINIASM ( input )
}

workflow test_dragonflye_skesa {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    DRAGONFLYE_RAVEN ( input )
}
