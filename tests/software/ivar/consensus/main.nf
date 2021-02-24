#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.save_mpileup = true
include { IVAR_CONSENSUS } from '../../../../software/ivar/consensus/main.nf' addParams( [ options: [args2: '-aa -A -d 0 -Q 0'] ] )

workflow test_ivar_consensus {
    def ref = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genomic.fasta", checkIfExists: true)
    def input = []
    input = [ [ id:'test'],
                file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.bam", checkIfExists: true) ]
    IVAR_CONSENSUS ( input, ref )
}
