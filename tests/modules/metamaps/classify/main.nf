#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR                } from '../../../../modules/untar/main.nf'
include { METAMAPS_MAPDIRECTLY } from '../../../../modules/metamaps/mapdirectly/main.nf'
include { METAMAPS_CLASSIFY    } from '../../../../modules/metamaps/classify/main.nf'

workflow test_metamaps_classify {

    classification_res = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    database_folder = [
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    UNTAR ( database )
    METAMAPS_MAPDIRECTLY ( input, UNTAR.out.untar )
    db_file = new File(UNTAR.out.untar, "/DB.fa")
    METAMAPS_CLASSIFY ( METAMAPS_MAPDIRECTLY, db_file )
}
