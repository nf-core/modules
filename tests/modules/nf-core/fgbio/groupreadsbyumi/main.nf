#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_GROUPREADSBYUMI } from '../../../../modules/fgbio/groupreadsbyumi/main.nf'

workflow test_fgbio_groupreadsbyumi {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_umi_unsorted_tagged_bam'], checkIfExists: true)
    ]

    FGBIO_GROUPREADSBYUMI ( input, 'Adjacency' )
}
