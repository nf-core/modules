#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = "./"

include { LONGRANGER_MKREF } from '../../../../modules/longranger/mkref/main.nf'


workflow test_longranger_mkref {
    LONGRANGER_MKREF([ [id : "pEimTen1"],
        file(params.tol_test_data['small_genome']['pEimTen1']['assembly']['canu_contigs_fasta'] , checkIfExists: true) ])
}
