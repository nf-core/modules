#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OPTITYPE } from '../../../modules/optitype/main.nf' addParams( options: ['args':'-e 1 -b 0.009', 'args2':'solver=glpk'] )

workflow test_optitype {
    input = [ [ id:'test', seq_type:'dna' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['example_hla_pe'], checkIfExists: true)
            ]

    OPTITYPE ( input )
}
