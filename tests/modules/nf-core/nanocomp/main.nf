#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOCOMP } from '../../../../modules/nf-core/nanocomp/main.nf'

workflow test_nanocomp_fastq {
    
    inputfile = Channel.fromPath(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    inputchannel = inputfile.toList()
    NANOCOMP ( "fastq", inputchannel )
}

workflow test_nanocomp_summary {
    
    inputfile = Channel.fromPath(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true)
    inputchannel = inputfile.toList()
    NANOCOMP ( "summary", inputchannel )
}
