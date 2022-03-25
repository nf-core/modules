#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LONGRANGER_MKREF } from '../../../../modules/longranger/mkref/main.nf'

workflow test_longranger_mkref {
    LONGRANGER_MKREF([ [id : "longranger test"],
        file('https://darwin.cog.sanger.ac.uk/longranger_nf_test/pEimTen1.contigs.fasta', checkIfExists: true) ])
}
