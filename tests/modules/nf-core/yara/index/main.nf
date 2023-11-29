#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YARA_INDEX } from '../../../../../modules/nf-core/yara/index/main.nf'

workflow test_yara_index {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    meta = [ id:'test' ]

    YARA_INDEX ( [ meta, input ] )
}
