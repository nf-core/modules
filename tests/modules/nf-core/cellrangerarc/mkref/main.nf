#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERARC_MKREF } from '../../../../../modules/nf-core/cellrangerarc/mkref/main.nf'
include { UNZIP } from '../../../../../modules/nf-core/unzip/main.nf'

workflow test_cellrangerarc_mkref {

    // fasta = file(params.test_data['mus_musculus']['genome']['genome_fasta'], checkIfExists: true)
    // gtf = file(params.test_data['mus_musculus']['genome']['genome_gtf'], checkIfExists: true)
    // motifs = file(params.test_data['mus_musculus']['genome']['genome_motif'], checkIfExists: true)
    // reference_config = file(params.test_data['mus_musculus']['genome']['genome_config'], checkIfExists: true)
    // reference_name = "cellrangerarc_reference"

    fasta = [[], file("/home/florian/Downloads/cellranger_testdata/chr19.fa.gz", checkIfExists: true)]
    gtf = file("/home/florian/Downloads/cellranger_testdata/chr19.filtered.gtf.gz", checkIfExists: true)
    motifs = file("/home/florian/Downloads/cellranger_testdata/genome_motifs.txt", checkIfExists: true)
    reference_config = file("/home/florian/Downloads/cellranger_testdata/genome_config.json", checkIfExists: true)
    reference_name = "cellrangerarc_reference"

    UNZIP( fasta )

    CELLRANGERARC_MKREF ( UNZIP.out.unzipped_archive.map { it[1] } + "/chr19.fa",
                            gtf,
                            motifs,
                            reference_config,
                            reference_name )
}
