#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MODKIT_PILEUP } from '../../../../../modules/nf-core/modkit/pileup/main.nf'

workflow test_modkit_pileup {
    
    bam_input = [
        [ id:'test', single_end:false ],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/nanopore/bam/bc_anchored_10_reads.sorted.bam', checkIfExists: true)
    ]
    bai_input = [
    file('https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/nanopore/bam/bc_anchored_10_reads.sorted.bam.bai', checkIfExists: true) ]
    output_bed = [ file("test.bed") ]

    MODKIT_PILEUP ( bam_input, bai_input, output_bed )
}