#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CLIPPY } from '../../../../modules/nf-core/clippy/main.nf'

workflow test_clippy {
    

    input = [
        [  id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_multi_interval_bed'], checkIfExists: true)
    ]

    CLIPPY ( 
        input,
        file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
        )
}
