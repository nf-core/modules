#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ_FILTER } from '../../../../software/lofreq/filter/main.nf' addParams( options: [:] )

workflow test_lofreq_filter {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]

    LOFREQ_FILTER ( input )
}
