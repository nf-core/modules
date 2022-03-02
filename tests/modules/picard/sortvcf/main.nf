#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_SORTVCF } from '../../../../modules/picard/sortvcf/main.nf'

workflow test_picard_sortvcf {

    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)
            ]

    PICARD_SORTVCF ( input )
}
