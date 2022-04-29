#!/usr/bin/env nextflow



include { LOFREQ_FILTER } from '../../../../modules/lofreq/filter/main.nf'

workflow test_lofreq_filter {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true) ]

    LOFREQ_FILTER ( input )
}
