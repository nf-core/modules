#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMOOTHXG_ITERATE } from '../../../../../modules/nf-core/smoothxg/iterate/main.nf'

workflow test_smoothxg_iterate {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['assembly_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG_ITERATE ( input )
}

workflow test_smoothxg_iterate_pangenomics {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['homo_sapiens']['pangenome']['pangenome_seqwish_gfa'], checkIfExists: true) ]
            ]

    SMOOTHXG_ITERATE ( input )
}
