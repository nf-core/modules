#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURGEDUPS_FASTK } from '../../../../modules/purgedups/fastk/main.nf'

workflow test_purgedups_fastk {
    
        data = [
        [ id:'test'], // meta map
        [
           file('/lustre/scratch123/tol/teams/grit/yy5/nf_test_small/ilDiaRubi1.trimmedReads.fasta.gz', checkIfExists: true),
        ],
    ]

    kmer=31
    myoutdir='/lustre/scratch123/tol/teams/grit/yy5/nf_test_small'
    FASTKDB='FASTKDB'

    PURGEDUPS_FASTK ( data,kmer,myoutdir,FASTKDB )
}
