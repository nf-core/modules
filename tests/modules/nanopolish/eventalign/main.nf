#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP                } from '../../../../modules/gunzip/main.nf'                addParams( options: [:] )
include { NANOPOLISH_EVENTALIGN } from '../../../../modules/nanopolish/eventalign/main.nf' addParams( options: [:] )
include { UNTAR                 } from '../../../../modules/untar/main.nf'                 addParams( options: [:] )

workflow test_nanopolish_eventalign {

    // input files
    input   = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['nanopore']['test_sorted_bam'], checkIfExists: true),
                file(params.test_data['sarscov2']['nanopore']['test_sorted_bam_bai'], checkIfExists: true)
                ]
    fast5   = file(params.test_data['sarscov2']['nanopore']['fast5_tar_gz'], checkIfExists: true)
    fastq   = file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf     = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
    threads = 1

    // calling necessary modules
    UNTAR ( fast5 )
    GUNZIP ( fastq )
    NANOPOLISH_EVENTALIGN (input, UNTAR.out.untar , GUNZIP.out.gunzip, fasta, gtf, threads)
}
