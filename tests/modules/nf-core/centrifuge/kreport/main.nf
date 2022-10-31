#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR                } from '../../../../../modules/nf-core/untar/main.nf'
include { CENTRIFUGE_CENTRIFUGE } from '../../../../../modules/nf-core/centrifuge/centrifuge/main.nf'
include { CENTRIFUGE_KREPORT    } from '../../../../../modules/nf-core/centrifuge/kreport/main.nf'

workflow test_centrifuge_kreport_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    =  [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz', checkIfExists: true) ]

    ch_db = UNTAR ( db )
    CENTRIFUGE_CENTRIFUGE ( input, ch_db.untar.map{ it[1] }, false, false, false )
    CENTRIFUGE_KREPORT ( CENTRIFUGE_CENTRIFUGE.out.results, ch_db.untar.map{ it[1] } )
}

workflow test_centrifuge_kreport_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
     db    =  [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz', checkIfExists: true) ]

    ch_db = UNTAR ( db )
    CENTRIFUGE_CENTRIFUGE ( input, ch_db.untar.map{ it[1] }, false, false, false )
    CENTRIFUGE_KREPORT ( CENTRIFUGE_CENTRIFUGE.out.results, ch_db.untar.map{ it[1] } )
}

