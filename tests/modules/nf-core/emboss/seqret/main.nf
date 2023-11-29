#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMBOSS_SEQRET } from '../../../../../modules/nf-core/emboss/seqret/main.nf'
include { GUNZIP } from '../../../../modules/nf-core/gunzip/main.nf'

workflow test_emboss_seqret_gb2embl {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true)
    ]

    GUNZIP ( input )
    EMBOSS_SEQRET ( GUNZIP.out.gunzip, 'embl' )
}

workflow test_emboss_seqret_gb2gff {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true)
    ]

    GUNZIP ( input )
    EMBOSS_SEQRET ( GUNZIP.out.gunzip, 'gff' )
}

workflow test_emboss_seqret_gb2pir {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true)
    ]

    GUNZIP ( input )
    EMBOSS_SEQRET ( GUNZIP.out.gunzip, 'pir' )
}

workflow test_emboss_seqret_gb2fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true)
    ]

    GUNZIP ( input )
    EMBOSS_SEQRET ( GUNZIP.out.gunzip, 'fasta' )
}
