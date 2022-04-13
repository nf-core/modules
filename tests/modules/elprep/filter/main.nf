#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ELPREP_FILTER } from '../../../../modules/elprep/filter/main.nf'

workflow test_elprep_filter {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    ELPREP_FILTER ( input, [], [], [],  [], [], [], [], [], [], [], [])
}
