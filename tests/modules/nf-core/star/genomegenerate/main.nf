#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { STAR_GENOMEGENERATE } from "$moduleDir/modules/nf-core/star/genomegenerate/main.nf"

workflow test_star_genomegenerate {
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    STAR_GENOMEGENERATE ( fasta, gtf )
}
