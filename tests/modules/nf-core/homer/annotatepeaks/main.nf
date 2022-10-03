#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_ANNOTATEPEAKS } from '../../../../modules/homer/annotatepeaks/main.nf'

workflow test_homer_annotatepeaks {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
            ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)
    HOMER_ANNOTATEPEAKS ( input, fasta, gtf )
}
