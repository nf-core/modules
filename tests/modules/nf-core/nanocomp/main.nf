#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOCOMP } from '../../../../modules/nf-core/nanocomp/main.nf'

workflow test_nanocomp {
    
    inputfile = Channel.fromPath(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    //inputfile = Channel.fromPath( '/scratch/wolkp/nfcore/test-datasets/data/genomics/sarscov2/nanopore/fastq/*' )
    inputchannel = inputfile.toList()
    NANOCOMP ( "fastq", inputchannel )
}
