#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOTUS_DOWNLOADDB                } from '../../../../../modules/nf-core/motus/downloaddb/main.nf'
include { MOTUS_PROFILE                   } from '../../../../../modules/nf-core/motus/profile/main.nf'
include { MOTUS_MERGE                     } from '../../../../../modules/nf-core/motus/merge/main.nf'
include { MOTUS_MERGE as MOTUS_MERGE_BIOM } from '../../../../../modules/nf-core/motus/merge/main.nf'

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

    MOTUS_DOWNLOADDB ( file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py') )
    MOTUS_PROFILE ( input, MOTUS_DOWNLOADDB.out.db )
    MOTUS_MERGE ( MOTUS_PROFILE.out.out.map{ [[id:"test"], it[1]] }.groupTuple(), MOTUS_DOWNLOADDB.out.db, MOTUS_PROFILE.out.versions.first() )

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

    MOTUS_DOWNLOADDB ( file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py') )
    MOTUS_PROFILE ( input, MOTUS_DOWNLOADDB.out.db )
    MOTUS_MERGE_BIOM ( MOTUS_PROFILE.out.out.map{ [[id:"test"], it[1]] }.groupTuple(), MOTUS_DOWNLOADDB.out.db, MOTUS_PROFILE.out.versions.first() )

}

