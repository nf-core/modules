#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { FASTK_FASTK          } from "$moduleDir/modules/nf-core/fastk/fastk/main.nf"
include { MERQURYFK_PLOIDYPLOT } from "$moduleDir/modules/nf-core/merquryfk/ploidyplot/main.nf"

workflow test_merquryfk_ploidyplot_png {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
    MERQURYFK_PLOIDYPLOT ( FASTK_FASTK.out.hist
        .join( FASTK_FASTK.out.ktab )
    )
}

workflow test_merquryfk_ploidyplot_pdf {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
    MERQURYFK_PLOIDYPLOT (
        FASTK_FASTK.out.hist
            .join( FASTK_FASTK.out.ktab )
    )
}
