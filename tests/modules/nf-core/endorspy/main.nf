#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENDORSPY } from '../../../../modules/nf-core/endorspy/main.nf'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT1 } from '../../../../modules/nf-core/samtools/flagstat/main.nf'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT2 } from '../../../../modules/nf-core/samtools/flagstat/main.nf'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT3 } from '../../../../modules/nf-core/samtools/flagstat/main.nf'
include { SAMTOOLS_VIEW } from '../../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_INDEX } from '../../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX2 } from '../../../../modules/nf-core/samtools/index/main.nf'

workflow test_endorspy_raw {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]


    SAMTOOLS_FLAGSTAT1 ( input )

    ch_input_flagstat = SAMTOOLS_FLAGSTAT1.out.flagstat.map{[it[0], it[1], [], []]}
    ENDORSPY ( ch_input_flagstat )
}

workflow test_endorspy_raw_qualityfilter {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]


    SAMTOOLS_FLAGSTAT1 ( input )
    SAMTOOLS_VIEW ( input, [], [], [] )
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.bam )
    input2 = SAMTOOLS_VIEW.out.bam
       .mix(SAMTOOLS_INDEX.out.bai)
       .groupTuple(by:0)
       .map{
            def meta = it[0]
            def bam = it[1][0]
            def bai = it[1][1]

            [meta, bam, bai]
       }
    SAMTOOLS_FLAGSTAT2 ( input2 )
    ch_input_flagstat = SAMTOOLS_FLAGSTAT1.out.flagstat.join(SAMTOOLS_FLAGSTAT2.out.flagstat).map{[it[0], it[1], it[2], []]}
    ENDORSPY ( ch_input_flagstat )
}

workflow test_endorspy_raw_qualityfilter_dedup {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]


    SAMTOOLS_FLAGSTAT1 ( input )
    SAMTOOLS_VIEW ( input, [], [], [] )
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.bam )
    input2 = SAMTOOLS_VIEW.out.bam
       .mix(SAMTOOLS_INDEX.out.bai)
       .groupTuple(by:0)
       .map{
            def meta = it[0]
            def bam = it[1][0]
            def bai = it[1][1]

            [meta, bam, bai]
       }
    SAMTOOLS_FLAGSTAT2 ( input2 )
    SAMTOOLS_FLAGSTAT3 ( input2 )

    ch_input_flagstat = SAMTOOLS_FLAGSTAT1.out.flagstat.join(SAMTOOLS_FLAGSTAT2.out.flagstat).join(SAMTOOLS_FLAGSTAT3.out.flagstat)
    ENDORSPY ( ch_input_flagstat )
}

workflow test_endorspy_raw_dedup {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]


    SAMTOOLS_FLAGSTAT1 ( input )
    SAMTOOLS_VIEW ( input, [], [], [] )
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.bam )
    input2 = SAMTOOLS_VIEW.out.bam
       .mix(SAMTOOLS_INDEX.out.bai)
       .groupTuple(by:0)
       .map{
            def meta = it[0]
            def bam = it[1][0]
            def bai = it[1][1]

            [meta, bam, bai]
       }
    SAMTOOLS_FLAGSTAT2 ( input2 )
    ch_input_flagstat = SAMTOOLS_FLAGSTAT1.out.flagstat.join(SAMTOOLS_FLAGSTAT2.out.flagstat).map{[it[0], it[1], [], it[2]]}
    ENDORSPY ( ch_input_flagstat )
}
workflow test_endorspy_qualityfilter_dedup {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]


    SAMTOOLS_FLAGSTAT1 ( input )
    SAMTOOLS_VIEW ( input, [], [], [] )
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.bam )
    input2 = SAMTOOLS_VIEW.out.bam
       .mix(SAMTOOLS_INDEX.out.bai)
       .groupTuple(by:0)
       .map{
            def meta = it[0]
            def bam = it[1][0]
            def bai = it[1][1]

            [meta, bam, bai]
       }
    SAMTOOLS_FLAGSTAT2 ( input2 )
    ch_input_flagstat = SAMTOOLS_FLAGSTAT1.out.flagstat.join(SAMTOOLS_FLAGSTAT2.out.flagstat).map{[it[0], [], it[1], it[2]]}
    ENDORSPY ( ch_input_flagstat )
}
