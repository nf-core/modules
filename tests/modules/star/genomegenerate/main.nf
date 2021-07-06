#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../software/star/genomegenerate/main.nf' addParams( options: [:] )

workflow test_star_genomegenerate {
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    STAR_GENOMEGENERATE ( fasta, gtf )
}
