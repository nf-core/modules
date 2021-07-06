#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../../../../software/samtools/faidx/main.nf' addParams( options: [:] )

workflow test_samtools_faidx {
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    SAMTOOLS_FAIDX ( fasta  )
}
