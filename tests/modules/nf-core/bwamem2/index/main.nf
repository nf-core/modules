#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../../../../modules/nf-core/bwamem2/index/main.nf'

workflow test_bwamem2_index {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWAMEM2_INDEX ( [ [id:'test'], fasta] )
}
