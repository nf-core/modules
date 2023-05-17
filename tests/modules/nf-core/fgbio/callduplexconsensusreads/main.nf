#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FGBIO_CALLDUPLEXCONSENSUSREADS } from '../../../../../modules/nf-core/fgbio/callduplexconsensusreads/main.nf'

workflow test_fgbio_callduplexconsensusreads {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_duplex_umi_grouped_bam'], checkIfExists: true)
    ]

    FGBIO_CALLDUPLEXCONSENSUSREADS ( input )
}
