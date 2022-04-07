#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR      } from '../../../modules/untar/main.nf'
include { CENTRIFUGE } from '../../../modules/centrifuge/main.nf'

workflow test_centrifuge_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    db    =  [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz', checkIfExists: true) ]
    db_name = "minigut_cf"
    save_unaligned = true
    save_aligned = false
    sam_format = false

    UNTAR ( db )
    CENTRIFUGE ( input, UNTAR.out.untar.map{ it[1] },db_name, save_unaligned, save_aligned, sam_format )

}

workflow test_centrifuge_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
     db    =  [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz', checkIfExists: true) ]
     db_name = "minigut_cf"
     save_unaligned = true
     save_aligned = false
     sam_format = false

    UNTAR ( db )
    CENTRIFUGE ( input, UNTAR.out.untar.map{ it[1] }, db_name, save_unaligned, save_aligned, sam_format )


}
