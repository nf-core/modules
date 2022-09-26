#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_DEPTH } from '../../../modules/samtools/depth/main.nf'
include { SEXDETERRMINE } from '../../../modules/sexdeterrmine/main.nf'

workflow test_sexdeterrmine {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test3_single_end_markduplicates_sorted_bam'], checkIfExists: true) ]

    SAMTOOLS_DEPTH ( input )
    SEXDETERRMINE ( SAMTOOLS_DEPTH.out.tsv, [] )
}
