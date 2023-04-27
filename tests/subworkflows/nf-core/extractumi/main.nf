#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXTRACTUMI } from '../../../../subworkflows/nf-core/extractumi/main.nf'

workflow test_extractumi {
    
    // channels enable parralle: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel
    fastq = [
    [[id:'gene1', single_end:false], [file('test_nucleo/fastq/seracare_0-5_R1_001ad.fastq.gz'), file('test_nucleo/fastq/seracare_0-5_R2_001ad.fastq.gz')]],
    [[id:'gene2', single_end:false], [file('test_nucleo/fastq/seracare_0-5_R1_001ae.fastq.gz'), file('test_nucleo/fastq/seracare_0-5_R2_001ae.fastq.gz')]]
    ]
    fastq = ch_fastq = Channel.fromList(fastq)
    EXTRACTUMI ( fastq )
}
