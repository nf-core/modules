#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { MOTUS_DOWNLOADDB } from "$moduleDir/modules/nf-core/motus/downloaddb/main.nf"
include { MOTUS_PROFILE    } from "$moduleDir/modules/nf-core/motus/profile/main.nf"

workflow test_motus_profile_single_end {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    MOTUS_DOWNLOADDB(file('https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py'))

    MOTUS_PROFILE ( input,  MOTUS_DOWNLOADDB.out.db )
}

