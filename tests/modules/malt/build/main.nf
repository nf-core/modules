#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MALT_BUILD } from '../../../../modules/malt/build/main.nf' addParams( options: [:] )

workflow test_malt_build {

    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    MALT_BUILD ( fastas, map )
}
