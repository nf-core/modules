#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRESTO_FILTERSEQ } from '../../../../../modules/nf-core/presto/filterseq/main.nf'
include { GUNZIP as GUNZIP_R1 } from '../../../../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_R2 } from '../../../../../modules/nf-core/gunzip/main.nf'


workflow test_presto_filterseq {

    input_R1 = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_airrseq_1_fastq_gz'], checkIfExists: true)
    ]
    input_R2 = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_airrseq_2_fastq_gz'], checkIfExists: true)
    ]

    GUNZIP_R1(input_R1)
    GUNZIP_R2(input_R2)

    ch_gunzip_R1 = GUNZIP_R1.out.gunzip
                        .map{ it -> [it[0].id, it[0], it[1] ]}
    ch_gunzip_R2 = GUNZIP_R2.out.gunzip
                                .map{ it -> [it[0].id, it[0], it[1] ]}
    ch_merged = ch_gunzip_R2.mix(ch_gunzip_R1)
                            .groupTuple()
                            .dump()
                            .map{ it -> [it[1][0], it[2][0], it[2][1]]}

    PRESTO_FILTERSEQ ( ch_merged )
}
