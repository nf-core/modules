#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.save_mpileup = true

include { IVAR_VARIANTS } from '../../../../software/ivar/variants/main.nf' addParams([:])

workflow test_ivar_variants_no_gff {
    def ref = file("${launchDir}/tests/data/fasta/sarscov2/MN908947.3.fa", checkIfExists: true)
    def input = []
    input = [ [ id:'test'], 
                file("${launchDir}/tests/data/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref )
}

params.gff = file("${launchDir}/tests/data/gff/sarscov2/MN908947.3.gff3", checkIfExists: true)

workflow test_ivar_variants_with_gff {
    
    def ref = file("${launchDir}/tests/data/fasta/sarscov2/MN908947.3.fa", checkIfExists: true)
    def input = []
    input = [ [ id:'test'], 
                file("${launchDir}/tests/data/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref )
}
