#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR                 } from '../../../../../modules/nf-core/untar/main.nf'
include { CENTRIFUGE_CENTRIFUGE } from '../../../../../modules/nf-core/centrifuge/centrifuge/main.nf'

workflow test_centrifuge_centrifuge_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    =  [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz', checkIfExists: true) ]
    save_unaligned = true
    save_aligned = false

    UNTAR ( db )
    CENTRIFUGE_CENTRIFUGE ( input, UNTAR.out.untar.map{ it[1] }, save_unaligned, save_aligned )

}

workflow test_centrifuge_centrifuge_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
     db    =  [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz', checkIfExists: true) ]
     save_unaligned = true
     save_aligned = false

    UNTAR ( db )
    CENTRIFUGE_CENTRIFUGE ( input, UNTAR.out.untar.map{ it[1] }, save_unaligned, save_aligned )


}
