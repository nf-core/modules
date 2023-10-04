#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GANON_BUILDCUSTOM } from '../../../../../modules/nf-core/ganon/buildcustom/main.nf'
include { GANON_CLASSIFY    } from '../../../../../modules/nf-core/ganon/classify/main.nf'

workflow test_ganon_classify {

    input_db = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]

    GANON_BUILDCUSTOM ( input_db, [], [] )
    GANON_CLASSIFY    ( input, GANON_BUILDCUSTOM.out.db.map{it[1]} )
}

workflow test_ganon_classify_pe {

    input_db = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true),
        ]
    ]

    GANON_BUILDCUSTOM ( input_db, [], [] )
    GANON_CLASSIFY    ( input, GANON_BUILDCUSTOM.out.db.map{it[1]} )
}




