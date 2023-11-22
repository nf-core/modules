#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVNATOR_CNVNATOR as CNVNATOR_RD        } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_HIST      } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_STAT      } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_PARTITION } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_CALL      } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'

workflow test_cnvnator {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    CNVNATOR_RD ( input, [[:],[]], [[:],[]], [[:],[]] )
    CNVNATOR_HIST ( [[:],[],[]], CNVNATOR_RD.out.root, [[:],[]], [[:],[]] )
    CNVNATOR_STAT ( [[:],[],[]], CNVNATOR_HIST.out.root, [[:],[]], [[:],[]] )
    CNVNATOR_PARTITION ( [[:],[],[]], CNVNATOR_STAT.out.root, [[:],[]], [[:],[]] )
    CNVNATOR_CALL ( [[:],[],[]], CNVNATOR_PARTITION.out.root, [[:],[]], [[:],[]] )
}
