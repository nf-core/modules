#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ZIP } from '../../../../modules/nf-core/zip/main.nf'

workflow test_zip {

    files_to_zip = [[], file(params.test_data['sarscov2']['genome']['genome_fasta'])]

    ZIP ( files_to_zip )
}

