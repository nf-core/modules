#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_FILTERCONSENSUSREADS } from '../../../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'

workflow test_fgbio_filterconsensusreads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_duplex_consensus_bam'], checkIfExists: true)
    ]

    fasta               = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    min_reads           = '3'
    min_baseq           = '45'
    max_base_error_rate = '0.2'

    FGBIO_FILTERCONSENSUSREADS ( input, fasta, min_reads, min_baseq, max_base_error_rate )
}
