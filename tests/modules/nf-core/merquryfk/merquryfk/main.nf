#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTK_FASTK                } from '../../../../../modules/nf-core/fastk/fastk/main.nf'
include { FASTK_FASTK as FASTK1      } from '../../../../../modules/nf-core/fastk/fastk/main.nf'
include { FASTK_FASTK as FASTK2      } from '../../../../../modules/nf-core/fastk/fastk/main.nf'
include { MERQURYFK_MERQURYFK } from '../../../../../modules/nf-core/merquryfk/merquryfk/main.nf'

workflow test_merquryfk_merquryfk_png {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    assembly = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    haplotigs = [
        [ id:'test', single_end:true ], []
    ]

    FASTK_FASTK ( input )
    MERQURYFK_MERQURYFK ( FASTK_FASTK.out.hist
        .join( FASTK_FASTK.out.ktab )
        .join( Channel.value( assembly ) )
        .join( Channel.value( haplotigs ) ),
        [], []

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

    haplotigs = [
        [ id:'test', single_end:true ], []
    ]

    FASTK_FASTK ( input )
    MERQURYFK_MERQURYFK ( FASTK_FASTK.out.hist
        .join( FASTK_FASTK.out.ktab )
        .join( Channel.value( assembly ) )
        .join( Channel.value( haplotigs ) ),
        [], []
    )
}

workflow test_merquryfk_merquryfk_trio {

    input1 = [
        [ id:'test1', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    input2 = [
        [ id:'test2', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
    ]
    input3 = [
        [ id:'test3', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true)
    ]
    assembly1 = [
        [ id:'test1', single_end:true ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    haplotigs = [
        [ id:'test', single_end:true ], []
    ]

    FASTK_FASTK ( input1 )
    FASTK1 ( input2 )
    FASTK2 ( input3 )
    MERQURYFK_MERQURYFK ( FASTK_FASTK.out.hist
        .join( FASTK_FASTK.out.ktab )
        .join( Channel.value( assembly1 ) )
        .join( Channel.value( haplotigs ) ),
        FASTK1.out.ktab,
        FASTK2.out.ktab
    )
}