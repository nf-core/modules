#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRA_PIPELINE } from '../../../../modules/ultra/pipeline/main.nf'
include { GUNZIP         } from '../../../../modules/gunzip/main.nf'
include { GFFREAD        } from '../../../../modules/gffread/main.nf'

workflow test_ultra_pipeline {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]
    GUNZIP ( input )

    gtf    = file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    GFFREAD ( gtf )

    ULTRA_PIPELINE ( GUNZIP.out.gunzip, genome, GFFREAD.out.gtf )
}
