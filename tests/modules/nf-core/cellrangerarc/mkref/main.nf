#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_ARC_MKREF } from '../../../../../modules/nf-core/cellrangerarc/mkref/main.nf'

workflow test_cellranger_arc_mkref {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motif'], checkIfExists: true)
    reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true)
    reference_name = "cellranger_arc_reference"

    CELLRANGER_ARC_MKREF ( fasta,
                            gtf,
                            motifs,
                            reference_config,
                            reference_name )
}
