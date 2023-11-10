#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTA_CLEAN_FCS } from '../../../../subworkflows/nf-core/fasta_clean_fcs/main.nf'

workflow test_fasta_clean_fcs {

    input = [
        [ id:'test', single_end:false, taxid:'9606' ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    database = [
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.gxi"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.gxs"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.taxa.tsv"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.meta.jsonl"),
        file("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.blast_div.tsv.gz")
    ]

    FASTA_CLEAN_FCS ( input, database )
}
