#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_INTERVALFILE } from '../../../../../modules/nf-core/purecn/intervalfile/main.nf'

workflow test_purecn_intervalfile {

    bed_input = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    sequence = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    genome = Channel.value("hg19")

    PURECN_INTERVALFILE ( bed_input, sequence, genome )
}
