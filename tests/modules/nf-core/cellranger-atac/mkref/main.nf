#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_ATAC_MKREF } from '../../../../../modules/nf-core/cellranger-atac/mkref/main.nf'

workflow test_cellranger_atac_mkref {

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) //TODO has to be provided!
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true) //TODO has to be provided!
    motifs = file(params.test_data['homo_sapiens']['genome']['genome_motif'], checkIfExists: true) //TODO has to be provided!
    reference_config = file(params.test_data['homo_sapiens']['genome']['genome_config'], checkIfExists: true) //TODO has to be provided!
    reference_name = "cellranger_atac_reference"

    CELLRANGER_ATAC_MKREF ( fasta,
                            gtf,
                            motifs,
                            reference_config,
                            reference_name )
}
