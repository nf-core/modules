#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_REFERENCE } from '../../../../../modules/nf-core/cnvkit/reference/main.nf'

workflow test_cnvkit_reference {

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    targets     = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    antitargets = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_antitarget_bed'], checkIfExists: true)

    CNVKIT_REFERENCE ( fasta, targets, antitargets )
}
