#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FIXMATE } from '../../../../modules/samtools/fixmate/main.nf'

workflow test_samtools_fixmate {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SAMTOOLS_FIXMATE ( input )

}
