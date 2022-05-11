#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = "./"

include { LONGRANGER_MKREF } from '../../../../modules/longranger/mkref/main.nf'


workflow test_longranger_mkref {
    LONGRANGER_MKREF([ [id : "pEimTen1"],
        file('https://tolit.cog.sanger.ac.uk/test-data/Eimeria_tenella/working/pEimTen1.canu.20190919/pEimTen1.contigs.fasta', checkIfExists: true) ])
}
