#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../../modules/nf-core/star/genomegenerate/main.nf'

workflow test_star_genomegenerate {
    fasta = [
        [ id:'test_fasta' ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    ]
    gtf = [
        [ id:'test_gtf' ], // meta map
        [ file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true) ]
    ]
    STAR_GENOMEGENERATE ( fasta, gtf )
}
