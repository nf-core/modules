#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQ_ALIGN } from '../../../../subworkflows/nf-core/fastq_align/main.nf'

workflow test_fastq_align {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    FASTQ_ALIGN ( input )
}
