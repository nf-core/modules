#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SUBSEQ } from '../../../../software/seqtk/subseq/main.nf' addParams( options: ['suffix':'.filtered'] )

workflow test_seqtk_subseq {
    
    sequences = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    filter_list = file(params.test_data['sarscov2']['genome']['test_bed_gz'], checkIfExists: true)

    SEQTK_SUBSEQ ( sequences, filter_list )
}
