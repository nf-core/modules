#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../../modules/nf-core/star/genomegenerate/main.nf'

workflow test_star_genomegenerate {
    fasta = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_fasta']).map { it -> [[id:it.Name], it] }.collect()
    gtf = Channel.fromPath(params.test_data['homo_sapiens']['genome']['genome_gtf']).map { it -> [[id:it.Name], it] }.collect()
    STAR_GENOMEGENERATE ( fasta, gtf )
}
