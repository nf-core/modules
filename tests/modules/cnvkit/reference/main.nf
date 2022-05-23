#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVKIT_REFERENCE } from '../../../../modules/cnvkit/reference/main.nf'

workflow test_cnvkit_reference {

    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    targets     = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    antitargets = file("/Users/susanne/Documents/repos/forks/modules/test_antitarget/output/cnvkit/test.antitarget.bed", checkIfExists: true)

    CNVKIT_REFERENCE ( fasta, targets, antitargets )
}
