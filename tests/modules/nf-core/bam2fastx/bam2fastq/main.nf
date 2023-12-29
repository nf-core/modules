#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM2FASTX_BAM2FASTQ } from '../../../../../modules/nf-core/bam2fastx/bam2fastq/main.nf'

workflow test_bam2fastx_bam2fastq {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['alz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['pacbio']['alzpbi'], checkIfExists: true)
    ]

    BAM2FASTX_BAM2FASTQ ( input )
}
