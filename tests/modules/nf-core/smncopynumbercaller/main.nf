#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMNCOPYNUMBERCALLER } from '../../../../modules/nf-core/smncopynumbercaller/main.nf' addParams( options: [args: ''] )

workflow test_smncopynumbercaller {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_rna_paired_end_sorted_chr6_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_rna_paired_end_sorted_chr6_bam_bai'], checkIfExists: true),
    ]

    SMNCOPYNUMBERCALLER ( input )
}

