#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNAPALIGNER_INDEX } from '../../../../../modules/nf-core/snapaligner/index/main.nf'

workflow test_snapaligner_index {
    fasta = [
        [id:"test"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        [],
        [],
        []
    ]
    SNAPALIGNER_INDEX ( fasta )
}
