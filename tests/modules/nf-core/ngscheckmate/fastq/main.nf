#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE_BUILD } from '../../../../../modules/nf-core/bowtie/build/main.nf'
include { BEDTOOLS_MAKEWINDOWS } from '../../../../../modules/nf-core/bedtools/makewindows/main.nf'
include { NGSCHECKMATE_PATTERNGENERATOR } from '../../../../../modules/nf-core/ngscheckmate/patterngenerator/main.nf'
include { NGSCHECKMATE_FASTQ } from '../../../../../modules/nf-core/ngscheckmate/fastq/main.nf'
include { GAWK as GAWK_BED } from '../../../../../modules/nf-core/gawk/main.nf'
include { GAWK as GAWK_FAI } from '../../../../../modules/nf-core/gawk/main.nf'
workflow test_ngscheckmate_fastq {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
    ]

    fasta    = [
        [ id: 'sarscov2' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]

    fasta_fai    = [
        [ id: 'sarscov2' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
        ]

    bowtie_fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    // create a bed file representing the entire genome from the fasta index
    GAWK_FAI(fasta_fai, [])

    // split bed file into individual base pairs
    BEDTOOLS_MAKEWINDOWS(GAWK_FAI.out.output)

    // add additional columns to represent the ref/alt alleles (just constant here)
    GAWK_BED(BEDTOOLS_MAKEWINDOWS.out.bed, [])

    // generate the .pt file for ngscheckmate to use
    BOWTIE_BUILD ( bowtie_fasta )
    NGSCHECKMATE_PATTERNGENERATOR ( GAWK_BED.out.output, fasta, BOWTIE_BUILD.out.index )

    NGSCHECKMATE_FASTQ ( input, NGSCHECKMATE_PATTERNGENERATOR.out.pt )
}

workflow test_ngscheckmate_fastq_paired {

    input = [
        [ id:'test_paired', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta    = [
        [ id: 'sarscov2' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        ]

    fasta_fai    = [
        [ id: 'sarscov2' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
        ]

    bowtie_fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BOWTIE_BUILD ( bowtie_fasta )

    // create a bed file representing the entire genome from the fasta index
    GAWK_FAI(fasta_fai, [])

    // split bed file into individual base pairs
    BEDTOOLS_MAKEWINDOWS(GAWK_FAI.out.output)

    // add additional columns to represent the ref/alt alleles (just constant here)
    GAWK_BED(BEDTOOLS_MAKEWINDOWS.out.bed, [])

    // generate the .pt file for ngscheckmate to use
    BOWTIE_BUILD ( bowtie_fasta )
    NGSCHECKMATE_PATTERNGENERATOR ( GAWK_BED.out.output, fasta, BOWTIE_BUILD.out.index )

    NGSCHECKMATE_FASTQ ( input, NGSCHECKMATE_PATTERNGENERATOR.out.pt )
}
