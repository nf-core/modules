#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP                      } from '../../../../modules/nf-core/gunzip/main.nf'
include { NANOPOLISH_INDEX_EVENTALIGN } from '../../../../modules/nf-core/nanopolish/main.nf'
include { UNTAR                       } from '../../../../modules/nf-core/untar/main.nf'

workflow test_nanopolish {

    // input files
    fast5   = Channel.of([[], file(params.test_data['sarscov2']['nanopore']['fast5_tar_gz'], checkIfExists: true)])
    fastq   = Channel.of([[], file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)])
    fasta   = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf     = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    // calling necessary modules
    UNTAR ( fast5 )
    fast5=UNTAR.out.untar.map{it->it[1]}
    GUNZIP ( fastq )
    fastq=GUNZIP.out.gunzip.map{it->it[1]}

    input   = Channel.of ([ [ id:'test' ], // meta map
                            fasta,
                            gtf,
                            file(params.test_data['sarscov2']['nanopore']['test_sorted_bam'], checkIfExists: true),
                            file(params.test_data['sarscov2']['nanopore']['test_sorted_bam_bai'], checkIfExists: true)
                          ])
    input
        .combine(fast5)
        .combine(fastq)
        .map{it -> [it[0],it[1],it[2],it[5],it[6],it[3],it[4]]}
        .set {input}

    NANOPOLISH_INDEX_EVENTALIGN ( input )
}
