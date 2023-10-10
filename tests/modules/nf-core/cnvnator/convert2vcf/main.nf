#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CNVNATOR_CNVNATOR as CNVNATOR_RD        } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_HIST      } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_STAT      } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_PARTITION } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_CALL      } from '../../../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CONVERT2VCF                    } from '../../../../../modules/nf-core/cnvnator/convert2vcf/main.nf'

workflow test_cnvnator_convert2vcf {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    CNVNATOR_RD ( input, [[:],[]], [[:],[]], [[:],[]] )
    CNVNATOR_HIST ( [[:],[],[]], CNVNATOR_RD.out.pytor, [[:],[]], [[:],[]] )
    CNVNATOR_STAT ( [[:],[],[]], CNVNATOR_HIST.out.pytor, [[:],[]], [[:],[]] )
    CNVNATOR_PARTITION ( [[:],[],[]], CNVNATOR_STAT.out.pytor, [[:],[]], [[:],[]] )
    CNVNATOR_CALL ( [[:],[],[]], CNVNATOR_PARTITION.out.pytor, [[:],[]], [[:],[]] )
    CNVNATOR_CONVERT2VCF (CNVNATOR_CALL.out.tab)
}
