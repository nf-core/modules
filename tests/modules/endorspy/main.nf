#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENDORSPY } from '../../../modules/endorspy/main.nf'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/samtools/flagstat/main.nf'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT2 } from '../../../modules/samtools/flagstat/main.nf'
include { SAMTOOLS_VIEW } from '../../../modules/samtools/view/main.nf'
include { SAMTOOLS_INDEX } from '../../../modules/samtools/index/main.nf'

workflow test_endorspy {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    

    SAMTOOLS_FLAGSTAT ( input )
    SAMTOOLS_VIEW ( input, [] )
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
    ch_input_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat.mix(SAMTOOLS_FLAGSTAT2.out.flagstat).groupTuple()
    ENDORSPY ( ch_input_flagstat )
}
