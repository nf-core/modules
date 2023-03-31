#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERATAC_MKREF } from '../../../../../modules/nf-core/cellrangeratac/mkref/main.nf'

workflow test_cellrangeratac_mkref {

    // fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    // gtf = file(params.test_data['homo_sapiens']['genome']['genome_scATAC_gtf'], checkIfExists: true)
    // motifs = file(params.test_data['homo_sapiens']['genome']['genome_motif'], checkIfExists: true)
    // reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true)
    fasta = file("/home/florian/Downloads/scatac_testdata/genome.fasta", checkIfExists: true)
    gtf = file("/home/florian/Downloads/scatac_testdata/genome.gtf", checkIfExists: true)
    motifs = file("/home/florian/Downloads/scatac_testdata/genome_motifs.txt", checkIfExists: true)
    reference_config = file("/home/florian/Downloads/scatac_testdata/config", checkIfExists: true)
    reference_name = "cellrangeratac_reference"

    CELLRANGERATAC_MKREF ( fasta,
                            gtf,
                            motifs,
                            reference_config,
                            reference_name )
}
