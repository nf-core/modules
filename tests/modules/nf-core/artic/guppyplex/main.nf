#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARTIC_GUPPYPLEX } from '../../../../modules/artic/guppyplex/main.nf'

process STAGE_FASTQ_DIR {
    input:
    tuple val(meta), path(fastq_file)

    output:
    tuple val(meta), path('fastq'), emit: fastq_dir

    script:
    """
    mkdir fastq
    mv ${fastq_file} fastq
    """
}

workflow test_artic_guppyplex {

    input =  [ [ id:'test' ],
               file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]

    STAGE_FASTQ_DIR ( input )

    ARTIC_GUPPYPLEX ( STAGE_FASTQ_DIR.out.fastq_dir )
}
