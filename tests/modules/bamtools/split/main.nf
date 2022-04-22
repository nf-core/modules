#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMTOOLS_SPLIT } from '../../../../modules/bamtools/split/main.nf'

workflow test_bamtools_split_single_input {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    BAMTOOLS_SPLIT ( input )
}

workflow test_bamtools_split_multiple_input {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)]
    [

    BAMTOOLS_SPLIT ( input )
}
