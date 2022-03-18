#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR        } from '../../../../modules/untar/main.nf'
include { ARTIC_MINION } from '../../../../modules/artic/minion/main.nf'

workflow test_artic_minion {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true)
    ]
    fast5_tar          = [ [], file(params.test_data['sarscov2']['nanopore']['fast5_tar_gz'], checkIfExists: true) ]
    sequencing_summary = file(params.test_data['sarscov2']['nanopore']['test_sequencing_summary'], checkIfExists: true)
    fasta              = file('https://github.com/artic-network/primer-schemes/raw/master/nCoV-2019/V3/nCoV-2019.reference.fasta', checkIfExists: true)
    bed                = file('https://github.com/artic-network/primer-schemes/raw/master/nCoV-2019/V3/nCoV-2019.primer.bed', checkIfExists: true)

    fast5_dir = UNTAR ( fast5_tar ).untar.map{ it[1] }

    ARTIC_MINION ( input, fast5_dir, sequencing_summary, fasta, bed, [], '', 'nCoV-2019', '3')
}
