#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FCS_FCSGX } from '../../../../../modules/nf-core/fcs/fcsgx/main.nf'

workflow test_fcs_fcsgx {

    input = [
        [ id:'test', taxid:'9606' ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true),
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
