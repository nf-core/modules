#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { GUNZIP } from "$moduleDir/modules/nf-core/gunzip/main.nf"
include { PRODIGAL } from "$moduleDir/modules/nf-core/prodigal/main.nf"
include { DEEPBGC_DOWNLOAD } from "$moduleDir/modules/nf-core/deepbgc/download/main.nf"
include { DEEPBGC_PIPELINE } from "$moduleDir/modules/nf-core/deepbgc/pipeline/main.nf"

workflow test_deepbgc_pipeline_gbk {

    input = [
        [ id:'test_gbk', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    DEEPBGC_DOWNLOAD ()
    GUNZIP ( input )
    PRODIGAL ( GUNZIP.out.gunzip, 'gbk' )
    DEEPBGC_PIPELINE ( PRODIGAL.out.gene_annotations, DEEPBGC_DOWNLOAD.out.db )
}

workflow test_deepbgc_pipeline_fa {

    input = [
        [ id:'test_fa', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    DEEPBGC_DOWNLOAD ()
    GUNZIP ( input )
    DEEPBGC_PIPELINE ( GUNZIP.out.gunzip, DEEPBGC_DOWNLOAD.out.db )
}
