#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FLYE } from '../../../modules/flye/main.nf' addParams( options: [:] )

workflow test_flye_pacbio_raw {
    input = [ [ id:'test', mode:'pacbio-raw' ], // meta map
              [ file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    FLYE ( input )
}

workflow test_flye_pacbio_corr {
    input = [ [ id:'test', mode:'pacbio-corr' ], // meta map
              [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
            ]

    FLYE ( input )
}

workflow test_flye_pacbio_hifi {
    input = [ [ id:'test', mode:'pacbio-hifi' ], // meta map
              [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
            ]

    FLYE ( input )
}

workflow test_flye_nano_raw {
    input = [ [ id:'test', mode:'nano-raw' ], // meta map
              [ file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    FLYE ( input )
}

workflow test_flye_nano_corr {
    input = [ [ id:'test', mode:'nano-corr' ], // meta map
              [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
            ]

    FLYE ( input )
}

workflow test_flye_nano_hq {
    input = [ [ id:'test', mode:'nano-hq' ], // meta map
              [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
            ]

    FLYE ( input )
}
