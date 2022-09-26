#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTK_FASTK as FASTK1 } from '../../../../../modules/nf-core/fastk/fastk/main.nf'
include { FASTK_FASTK as FASTK2 } from '../../../../../modules/nf-core/fastk/fastk/main.nf'
include { MERQURYFK_KATCOMP } from '../../../../../modules/nf-core/merquryfk/katcomp/main.nf'

workflow test_merquryfk_katcomp_png {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK1 ( input )
    FASTK2 ( input )
    MERQURYFK_KATCOMP (
        FASTK1.out.hist
            .join( FASTK1.out.ktab )
            .join( FASTK2.out.hist )
            .join( FASTK2.out.ktab )
    )
}

workflow test_merquryfk_katcomp_pdf {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK1 ( input )
    FASTK2 ( input )
    MERQURYFK_KATCOMP (
        FASTK1.out.hist
            .join( FASTK1.out.ktab )
            .join( FASTK2.out.hist )
            .join( FASTK2.out.ktab )
    )
}
