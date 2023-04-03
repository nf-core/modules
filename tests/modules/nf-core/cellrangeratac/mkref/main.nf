#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGERATAC_MKREF } from '../../../../../modules/nf-core/cellrangeratac/mkref/main.nf'
include { UNZIP } from '../../../../../modules/nf-core/unzip/main.nf'

workflow test_cellrangeratac_mkref {

    fasta = [ [], file(params.test_data['homo_sapiens']['genome']['genome_1_fasta'], checkIfExists: true) ]
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_1_gtf'], checkIfExists: true)
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motifs'], checkIfExists: true)
    reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true)
    reference_name = "cellrangeratac_reference"

    UNZIP( fasta )

    CELLRANGERATAC_MKREF ( UNZIP.out.unzipped_archive.map { it[1] } + "/genome.fasta",
                            gtf,
                            motifs,
                            reference_config,
                            reference_name )
}
