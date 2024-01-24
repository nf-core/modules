#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PBTK_BAM2FASTQ } from '../../../../../modules/nf-core/pbtk/bam2fastq/main.nf'

workflow test_pbtk_bam2fastq {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['pacbio']['alz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['pacbio']['alzpbi'], checkIfExists: true),
    ]

    PBTK_BAM2FASTQ ( input )
}
