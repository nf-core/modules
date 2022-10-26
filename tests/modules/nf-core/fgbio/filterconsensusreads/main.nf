#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_FILTERCONSENSUSREADS } from '../../../../../modules/nf-core/fgbio/filterconsensusreads/main.nf'

workflow test_fgbio_filterconsensusreads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_duplex_consensus_bam'], checkIfExists: true)
    ]

    fasta               = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    FGBIO_FILTERCONSENSUSREADS ( input, fasta )
}
