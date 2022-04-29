#!/usr/bin/env nextflow



include { BLAST_MAKEBLASTDB } from '../../../../modules/blast/makeblastdb/main.nf'
include { BLAST_BLASTN      } from '../../../../modules/blast/blastn/main.nf'

workflow test_blast_blastn {
    input = [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( input )
    BLAST_BLASTN ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db )
}
