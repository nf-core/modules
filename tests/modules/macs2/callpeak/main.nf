#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MACS2_CALLPEAK } from '../../../../modules/macs2/callpeak/main.nf'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_CTRL    } from '../../../../modules/macs2/callpeak/main.nf'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_BED     } from '../../../../modules/macs2/callpeak/main.nf'

workflow test_macs2_callpeak_bed {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['pacbio']['genemodel1'], checkIfExists: true)],
                  []]

    MACS2_CALLPEAK_BED ( input, 4000 )
}

workflow test_macs2_callpeak {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  []]

    MACS2_CALLPEAK ( input, 40000 )
}

workflow test_macs2_callpeak_ctrl {
    input     = [ [ id:'test', single_end:false ], // meta map
                  [ file( params.test_data['homo_sapiens']['illumina']['test_paired_end_name_sorted_bam'], checkIfExists: true) ],
                  [ file( params.test_data['homo_sapiens']['illumina']['test2_paired_end_name_sorted_bam'], checkIfExists: true) ]]

    MACS2_CALLPEAK_CTRL ( input, 40000 )
}
