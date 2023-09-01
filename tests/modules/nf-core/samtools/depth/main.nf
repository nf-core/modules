#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_DEPTH } from '../../../../../modules/nf-core/samtools/depth/main.nf'

workflow test_samtools_depth {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true) ]

    intervals = [ [ id:'bed' ],
                  file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true) ]

    SAMTOOLS_DEPTH ( input, intervals )
}
