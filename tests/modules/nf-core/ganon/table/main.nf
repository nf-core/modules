#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GANON_BUILDCUSTOM } from '../../../../../modules/nf-core/ganon/buildcustom/main.nf'
include { GANON_CLASSIFY    } from '../../../../../modules/nf-core/ganon/classify/main.nf'
include { GANON_REPORT      } from '../../../../../modules/nf-core/ganon/report/main.nf'
include { GANON_TABLE       } from '../../../../../modules/nf-core/ganon/table/main.nf'

workflow test_ganon_classify {

    input_db = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    input = Channel.fromList([
        [ [ id:'test', single_end:true ], file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ],
        [ [ id:'test2', single_end:false ], [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true), file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ] ]
    ])

    GANON_BUILDCUSTOM ( input_db, [], []                                               )
    GANON_CLASSIFY    ( input                    , GANON_BUILDCUSTOM.out.db.map{it[1]} )
    GANON_REPORT      ( GANON_CLASSIFY.out.report, GANON_BUILDCUSTOM.out.db.map{it[1]} )
    GANON_TABLE       ( GANON_REPORT.out.tre.collect{it[1]}.map{[[id: "db1"], it]}        )

}
