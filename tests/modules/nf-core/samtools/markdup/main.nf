#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_COLLATE } from '../../../../../modules/nf-core/samtools/collate/main.nf'
include { SAMTOOLS_FIXMATE } from '../../../../../modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_SORT    } from '../../../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_MARKDUP } from '../../../../../modules/nf-core/samtools/markdup/main.nf'

workflow test_samtools_markdup {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SAMTOOLS_COLLATE ( input, [] )
    SAMTOOLS_FIXMATE ( SAMTOOLS_COLLATE.out.bam )
    SAMTOOLS_SORT ( SAMTOOLS_FIXMATE.out.bam )
    SAMTOOLS_MARKDUP ( SAMTOOLS_SORT.out.bam, [] )

}
