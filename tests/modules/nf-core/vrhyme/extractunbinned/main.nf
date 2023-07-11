#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP                    } from '../../../../../modules/nf-core/gunzip/main.nf'
include { VRHYME_VRHYME             } from '../../../../../modules/nf-core/vrhyme/vrhyme/main.nf'
include { VRHYME_EXTRACTUNBINNED    } from '../../../../../modules/nf-core/vrhyme/extractunbinned/main.nf'

workflow test_vrhyme_extractunbinned {

    reads = [
        [ id:'test', single_end:false ], // meta map
        [
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    GUNZIP ( fasta )

    VRHYME_VRHYME ( reads , GUNZIP.out.gunzip )

    VRHYME_EXTRACTUNBINNED ( VRHYME_VRHYME.out.membership , GUNZIP.out.gunzip )
}
