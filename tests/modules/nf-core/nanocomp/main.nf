#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NANOCOMP } from '../../../../modules/nf-core/nanocomp/main.nf'

workflow test_nanocomp_fastq {
    
    input = [[id: "test_"],[file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)]]
    NANOCOMP ( input )
}

workflow test_nanocomp_summary {
    input = [[id: "test_"],[file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true)]]
    NANOCOMP ( input )
}