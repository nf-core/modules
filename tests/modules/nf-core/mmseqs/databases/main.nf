#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MMSEQS_DATABASES } from '../../../../../modules/nf-core/mmseqs/databases/main.nf'

workflow test_mmseqs_databases {

    input = "SILVA"

    MMSEQS_DATABASES ( input )
}
