#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNAPALIGNER_INDEX } from '../../../../modules/snapaligner/index/main.nf'

workflow test_snapaligner_index {
    SNAPALIGNER_INDEX ( file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),[],[],[])
}
