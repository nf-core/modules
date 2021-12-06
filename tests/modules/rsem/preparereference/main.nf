#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RSEM_PREPAREREFERENCE }   from '../../../../modules/rsem/preparereference/main.nf'

workflow test_rsem_preparereference {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    RSEM_PREPAREREFERENCE ( fasta, gtf )
}
