#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FCS_FCSGX } from '../../../../../modules/nf-core/fcs/fcsgx/main.nf'

workflow test_fcs_fcsgx {

    input = [
        [ id:'test', taxid:'6973' ], // meta map
        file("/Users/tillenglert/Documents/FCSGX/fcsgx_test.fa.gz"), // file(params.test_data['Blattella_germanica']['genome']['genome_fa_gz'], checkIfExists: true),
    ]

    database = [
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.gxi"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.gxs"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.taxa.tsv"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.meta.jsonl"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.blast_div.tsv.gz")
    ]
    FCS_FCSGX ( input, database )
}
