#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ELPREP_SPLIT } from '../../../../../modules/nf-core/elprep/split/main.nf'
include { ELPREP_MERGE } from '../../../../../modules/nf-core/elprep/merge/main.nf'

workflow test_elprep_merge {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    ELPREP_SPLIT ( input )
    ELPREP_MERGE ( ELPREP_SPLIT.out.bam )
}
