#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YARA_INDEX } from '../../../../modules/yara/index/main.nf'

workflow test_yara_index {

    input = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    YARA_INDEX ( input )
}
