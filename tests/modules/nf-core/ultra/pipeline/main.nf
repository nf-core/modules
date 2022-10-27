#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { ULTRA_PIPELINE } from "$moduleDir/modules/nf-core/ultra/pipeline/main.nf"
include { GUNZIP         } from "$moduleDir/modules/nf-core/gunzip/main.nf"
include { GFFREAD        } from "$moduleDir/modules/nf-core/gffread/main.nf"

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
