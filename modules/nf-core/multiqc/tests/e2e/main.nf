#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from '../../../fastqc/main'
include { BOWTIE_ALIGN } from '../../../bowtie/align/main'
include { MULTIQC } from "../../main"

workflow test_multiqc {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    FASTQC  ( input )
    BOWTIE_BUILD ( fasta )
    BOWTIE_ALIGN ( input, BOWTIE_BUILD.out.index )

    ch_multiqc_files = Channel.empty()
        .mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        .mix(BOWTIE_ALIGN.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC ( ch_multiqc_files.collect(), [], [], [] )
}
