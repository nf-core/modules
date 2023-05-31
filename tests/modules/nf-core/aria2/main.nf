#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARIA2 } from '../../../../modules/nf-core/aria2/main.nf'

workflow test_aria2 {

    input     = params.test_data['sarscov2']['illumina']['test_single_end_bam']

    ARIA2 ( input )
}
