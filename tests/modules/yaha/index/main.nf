#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { YAHA_INDEX } from '../../../../modules/yaha/index/main.nf'

workflow test_yaha_index {

    fasta = Channel.of([[id: 'homo_sapien'], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)])

    YAHA_INDEX ( fasta )
}
