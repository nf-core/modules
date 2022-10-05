#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMTNOTE } from '../../../../modules/nf-core/hmtnote/main.nf'

workflow test_hmtnote {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
            ]

    HMTNOTE ( input)
}
