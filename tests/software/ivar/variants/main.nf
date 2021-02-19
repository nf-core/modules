#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_VARIANTS } from '../../../../software/ivar/variants/main.nf' addParams([:])

workflow test_ivar_variants_no_gff_no_mpileup {
    params.gff          = false
    params.save_mpileup = false

    def ref   = file("${launchDir}/tests/data/fasta/sarscov2/MN908947.3.fa", checkIfExists: true)
    def dummy = file("${launchDir}/tests/data/dummy/dummy_file.txt", checkIfExists: true)
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref, dummy )
}

workflow test_ivar_variants_no_gff_with_mpileup {
    params.gff          = false
    params.save_mpileup = true

    def ref   = file("${launchDir}/tests/data/fasta/sarscov2/MN908947.3.fa", checkIfExists: true)
    def dummy = file("${launchDir}/tests/data/dummy/dummy_file.txt", checkIfExists: true)
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref, dummy )
}

workflow test_ivar_variants_with_gff_with_mpileup {
    params.gff          = true
    params.save_mpileup = true

    def ref = file("${launchDir}/tests/data/fasta/sarscov2/MN908947.3.fa", checkIfExists: true)
    def gff = file("${launchDir}/tests/data/gff/sarscov2/MN908947.3.gff3", checkIfExists: true)
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref, gff )
}
