#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMOOTHXG_ONCE } from '../../../../../modules/nf-core/smoothxg/once/main.nf'

workflow test_smoothxg_once {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG_ONCE ( input )
}

workflow test_smoothxg_once_pangenomics {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_seqwish_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG_ONCE ( input )
}
