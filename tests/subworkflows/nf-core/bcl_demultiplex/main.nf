#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCL_DEMULTIPLEX } from '../../../../subworkflows/nf-core/bcl_demultiplex/main.nf'

workflow test_bcl_demultiplex_bclconvert {

    ch_flowcell = Channel.value([
            [id:'test', lane:1 ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell_samplesheet'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)])

    BCL_DEMULTIPLEX ( ch_flowcell, "bclconvert" )
}

workflow test_bcl_demultiplex_bcl2fastq {

    ch_flowcell = Channel.value([
            [id:'test', lane:1 ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell_samplesheet'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_flowcell'], checkIfExists: true)])

    BCL_DEMULTIPLEX ( ch_flowcell, "bcl2fastq" )
}
