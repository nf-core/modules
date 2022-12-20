#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_GTFTOGENEPRED          } from '../../../../../modules/nf-core/ucsc/gtftogenepred/main.nf'
include { PICARD_COLLECTRNASEQMETRICS } from '../../../../../modules/nf-core/picard/collectrnaseqmetrics/main.nf'

workflow test_picard_collectrnaseqmetrics {

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gtf   = file(params.test_data['sarscov2']['genome']['genome_gtf'], checkIfExists: true)

    input_gtftogenepred = [
        [ id:'test'],
        gtf
    ]

    UCSC_GTFTOGENEPRED(input_gtftogenepred)
    
    input_collectrnaseqmetrics = [
        [ id:'test', single_end:false, strandedness:'forward' ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PICARD_COLLECTRNASEQMETRICS(
        input_collectrnaseqmetrics,
        UCSC_GTFTOGENEPRED.out.refflat.map{ it[1] },
        fasta,
        []
    )

}
