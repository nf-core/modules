#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRA_PIPELINE } from '../../../../modules/ultra/pipeline/main.nf' addParams( options: [:] )
include { GUNZIP }         from '../../../../modules/gunzip/main.nf'         addParams( options: [:] )
include { GFFREAD }        from '../../../../modules/gffread/main.nf'        addParams( options: [suffix: "_sorted"] )

workflow test_ultra_pipeline {

    fastq  = file(params.test_data['homo_sapiens']['pacbio']['hifi']        , checkIfExists: true)
    gtf    = file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    GUNZIP(fastq)
    GFFREAD(gtf)

    ULTRA_PIPELINE (
        [ [ id:'test', single_end:false ], GUNZIP.out.gunzip ],
        genome,
        GFFREAD.out.gtf )
}
