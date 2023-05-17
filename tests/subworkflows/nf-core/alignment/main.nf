#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGNMENT } from '../../../../subworkflows/nf-core/alignment/main.nf'

workflow test_alignment {
    
    // channels enable parralle: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel
    fastqs = [
    [[id:'gene', single_end:false], [file('output/fastp/merged_fastq_1.fastp.fastq.gz'), file('output/fastp/merged_fastq_2.fastp.fastq.gz')]]
    ]
    reference = [
        [id:'reference'], 
        file('test_nucleo/reference/')
    ]
    fastqs = ch_fastq = Channel.fromList(fastqs)
    ALIGNMENT ( fastqs, reference)
}
