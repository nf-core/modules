#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMTNOTE_ANNOTATE } from '../../../../modules/nf-core/hmtnote/annotate/main.nf'

workflow test_hmtnote_annotate {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
            ]

    HMTNOTE_ANNOTATE ( input)
}
