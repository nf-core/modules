#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MOTUS_DOWNLOADDB } from '../../../../modules/motus/downloaddb/main.nf'
include { MOTUS_PROFILE    } from '../../../../modules/motus/profile/main.nf'

workflow test_motus_profile_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MOTUS_DOWNLOADDB(file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py'))

    MOTUS_PROFILE ( input,  MOTUS_DOWNLOADDB.out.db )
}

workflow test_motus_profile_paired_end {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
          file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
    ]

    MOTUS_DOWNLOADDB(file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py'))

    MOTUS_PROFILE ( input,  MOTUS_DOWNLOADDB.out.db)
}

workflow test_motus_profile_bam {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]

    MOTUS_DOWNLOADDB(file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py'))

    MOTUS_PROFILE ( input,  MOTUS_DOWNLOADDB.out.db)
}
