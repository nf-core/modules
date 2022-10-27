#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { FASTK_FASTK         } from "$moduleDir/modules/nf-core/fastk/fastk/main.nf"
include { MERQURYFK_MERQURYFK } from "$moduleDir/modules/nf-core/merquryfk/merquryfk/main.nf"

workflow test_merquryfk_merquryfk_png {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    assembly = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
    MERQURYFK_MERQURYFK ( FASTK_FASTK.out.hist
        .join( FASTK_FASTK.out.ktab )
        .join( Channel.value( assembly ) )
    )
}

workflow test_merquryfk_merquryfk_pdf {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    assembly = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    FASTK_FASTK ( input )
    MERQURYFK_MERQURYFK ( FASTK_FASTK.out.hist
        .join( FASTK_FASTK.out.ktab )
        .join( Channel.value( assembly ) )
    )
}
