#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_INTERVALFILE } from '../../../../../modules/nf-core/purecn/intervalfile/main.nf'

workflow test_purecn_intervalfile {

    meta = [id: "test"]
    bed_input = [meta, file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)]
    sequence = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    meta2 = [id: "fasta"]
    sequence_input = [meta2, sequence]
    genome = Channel.value("hg38")

    PURECN_INTERVALFILE ( bed_input, sequence_input, genome )
}
