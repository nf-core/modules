#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../software/blast/makeblastdb/main.nf' addParams( options: ['args': '-dbtype nucl'] )

workflow test_blast_makeblastdb {
    input =  [ file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
}
