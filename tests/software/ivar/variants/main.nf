#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IVAR_VARIANTS } from '../../../../software/ivar/variants/main.nf' addParams([:])

workflow test_ivar_variants_no_gff_no_mpileup {
    params.gff          = false
    params.save_mpileup = false

    def ref   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def dummy = file("dummy_file.txt")
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref, dummy )
}

workflow test_ivar_variants_no_gff_with_mpileup {
    params.gff          = false
    params.save_mpileup = true

    def ref   = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def dummy = file("dummy_file.txt")
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref, dummy )
}

workflow test_ivar_variants_with_gff_with_mpileup {
    params.gff          = true
    params.save_mpileup = true

    def ref = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def gff = file("${launchDir}/tests/data/genomics/sarscov2/gtf/test_genomic.gff3", checkIfExists: true)
    def input = []
    input = [ [ id:'test'],
              file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
    IVAR_VARIANTS ( input, ref, gff )
}
