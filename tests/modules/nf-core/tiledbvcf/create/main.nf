#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TILEDBVCF_CREATE } from '../../../../../modules/nf-core/tiledbvcf/create/main.nf'

workflow test_tiledbvcf_create {

    input = Channel.value('my_vcf_dataset')

    TILEDBVCF_CREATE ( input )
}
