#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.recursion = true

include { SMOOTHXG } from '../../../../modules/nf-core/smoothxg/main.nf'

workflow test_smoothxg {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG ( input )
}

workflow test_smoothxg_recurse_times {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG.recurse ( input ).times ( 2 )
}

workflow test_smoothxg_pangenomics {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_seqwish_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG ( input )
}
