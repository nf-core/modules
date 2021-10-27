#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../../modules/gunzip/main.nf' addParams( options: [:] )
include { NANOPOLISH_EVENTALIGN } from '../../../../modules/nanopolish/eventalign/main.nf' addParams( options: [:] )
include { UNTAR } from '../../../../modules/untar/main.nf' addParams( options: [:] )

workflow test_nanopolish_eventalign {
    input = [ [ id:'test' ], // meta map
                file(params.test_data['sarscov2']['nanopore']['test_sorted_bam'], checkIfExists: true)
                file(params.test_data['sarscov2']['nanopore']['test_sorted_bam_bai'], checkIfExists: true)
            ]

    fast5 = file(params.test_data['sarscov2']['nanopore']['fast5_tar_gz'], checkIfExists: true)
    UNTAR ( fast5 )

    fastq = file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    GUNZIP ( )

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    NANOPOLISH_EVENTALIGN ( input, UNTAR.out.untar , fasta, params.threads)
}
