#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX          } from '../../../../modules/nf-core/bwa/index/main.nf'
include { FASTQ_ALIGN_BWAALN } from '../../../../subworkflows/nf-core/fastq_align_bwaaln/main.nf'


workflow test_fastq_align_bwaaln_singleend {

    input = Channel.of([
            [ id:'test', single_end:true ], // meta map
            [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
        ])
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    FASTQ_ALIGN_BWAALN ( input, BWA_INDEX.out.index )
}

workflow test_fastq_align_bwaaln_paired_end {

    input = Channel.of([
                [ id:'test', single_end:false ], // meta map
                [
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                    file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                ]
    ])
            fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    FASTQ_ALIGN_BWAALN ( input, BWA_INDEX.out.index )

}

workflow test_fastq_align_bwaaln_both {

    input = Channel.fromList(
        [
            [ [ id:'test',  single_end:false ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ] ],
            [ [ id:'test2', single_end:true ],  [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ] ]
        ]
    )
    fasta = [
        [id: 'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BWA_INDEX ( fasta )
    FASTQ_ALIGN_BWAALN ( input, BWA_INDEX.out.index )

}

workflow test_fastq_align_bwa_multiref {

    input = Channel.fromList(
        [
            [ [ id:'test',  single_end:false ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ] ],
            [ [ id:'test2', single_end:true ],  [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ] ]
        ]
    )

    fasta = Channel.fromList( [
        [ [id: 'test'], file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ],
        [ [id: 'test2'], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    ] )

    BWA_INDEX ( fasta )
    FASTQ_ALIGN_BWAALN ( input, BWA_INDEX.out.index )
}
