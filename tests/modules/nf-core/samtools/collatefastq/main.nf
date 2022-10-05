#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_COLLATEFASTQ } from '../../../../../modules/nf-core/samtools/collatefastq/main.nf'

workflow test_samtools_collatefastq {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    SAMTOOLS_COLLATEFASTQ ( input )
}
