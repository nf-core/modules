#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOTUS_DOWNLOADDB                } from '../../../../modules/motus/downloaddb/main.nf'
include { MOTUS_PROFILE                   } from '../../../../modules/motus/profile/main.nf'
include { MOTUS_MERGE                     } from '../../../../modules/motus/merge/main.nf'
include { MOTUS_MERGE as MOTUS_MERGE_BIOM } from '../../../../modules/motus/merge/main.nf'

workflow test_motus_merge {

    input = Channel.fromList(
        [
            [
                [ id:'test', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ],
            [
                [ id:'test2', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]
        ]
    )

    //MOTUS_DOWNLOADDB ( file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py') )
    MOTUS_PROFILE ( input, file("/home/jfellows/.conda/envs/motus/lib/python3.9/site-packages/motus/db_mOTU", checkIfExists: true) )
    MOTUS_MERGE ( MOTUS_PROFILE.out.out.map{ [[id:"test"], it[1]] }.groupTuple(), file("/home/jfellows/.conda/envs/motus/lib/python3.9/site-packages/motus/db_mOTU", checkIfExists: true), MOTUS_PROFILE.out.versions.first(), false )

}

workflow test_motus_merge_biom {

    input = Channel.fromList(
        [
            [
                [ id:'test', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ],
            [
                [ id:'test2', single_end:true ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
            ]
        ]
    )

    //MOTUS_DOWNLOADDB ( file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py') )
    MOTUS_PROFILE ( input, file("/home/jfellows/.conda/envs/motus/lib/python3.9/site-packages/motus/db_mOTU", checkIfExists: true) )
    MOTUS_MERGE_BIOM ( MOTUS_PROFILE.out.out.map{ [[id:"test"], it[1]] }.groupTuple(), file("/home/jfellows/.conda/envs/motus/lib/python3.9/site-packages/motus/db_mOTU", checkIfExists: true), MOTUS_PROFILE.out.versions.first(), false )

}

