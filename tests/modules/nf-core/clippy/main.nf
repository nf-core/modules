#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLIPPY } from '../../../../modules/nf-core/clippy/main.nf'

workflow test_clippy {
    

    input = [
        [  id:'test' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/clipseq/crosslinks/clippy.bed", checkIfExists: true)
    ]

    CLIPPY ( 
        input,
        file(params.test_data['homo_sapiens']['genome']['genome_21_gencode_gtf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
        )
}
