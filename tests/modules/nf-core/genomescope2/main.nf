#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MERYL_COUNT     } from '../../../../modules/nf-core/meryl/count/main.nf'
include { MERYL_HISTOGRAM } from '../../../../modules/nf-core/meryl/histogram/main.nf'
include { GENOMESCOPE2    } from '../../../../modules/nf-core/genomescope2/main.nf'

workflow test_genomescope2 {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true)
    ]

    MERYL_COUNT ( input )
    MERYL_HISTOGRAM ( MERYL_COUNT.out.meryl_db )
    GENOMESCOPE2 ( MERYL_HISTOGRAM.out.hist )
}
