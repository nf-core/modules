#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DSHBIO_EXPORTSEGMENTS } from '../../../../modules/dshbio/exportsegments/main.nf'

workflow test_dshbio_exportsegments {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    DSHBIO_EXPORTSEGMENTS ( input )
}

workflow test_dshbio_exportsegments_bgz {
    input = [ [ id:'test_bgz' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa_bgz'], checkIfExists: true) ]
            ]

    DSHBIO_EXPORTSEGMENTS ( input )
}

workflow test_dshbio_exportsegments_gz {
    input = [ [ id:'test_gz' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa_gz'], checkIfExists: true) ]
            ]

    DSHBIO_EXPORTSEGMENTS ( input )
}

workflow test_dshbio_exportsegments_zst {
    input = [ [ id:'test_zst' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa_zst'], checkIfExists: true) ]
            ]

    DSHBIO_EXPORTSEGMENTS ( input )
}
