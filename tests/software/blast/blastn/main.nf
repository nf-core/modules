#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../software/blast/makeblastdb/main.nf' addParams( options: ['args': '-dbtype nucl'] )
include { BLAST_BLASTN } from '../../../../software/blast/blastn/main.nf' addParams( options: [:] )

workflow test_blast_blastn {

    def input = []
    input = [ file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true) ]

    BLAST_MAKEBLASTDB (input)
    BLAST_BLASTN ([ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db)
}
