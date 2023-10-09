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

    input1 = [
        [ id:'test', single_end:false ],
        file("/home/ramprasad.neethiraj/nextflow/raredisease/work/91/0e416924ada7507395408109876bcb/earlycasualcaiman_T1.bam", checkIfExists: true),
        file("/home/ramprasad.neethiraj/nextflow/raredisease/work/2a/b8547d48fb0115b3268f1d89dd66d2/earlycasualcaiman_T1.bam.bai", checkIfExists: true)
	]

    CNVNATOR_RD ( input1, [[:],[]], [], [] )
    CNVNATOR_HIST ( [[:],[],[]], CNVNATOR_RD.out.pytor, [], [] )
    CNVNATOR_STAT ( [[:],[],[]], CNVNATOR_HIST.out.pytor, [], [] )
    CNVNATOR_PARTITION ( [[:],[],[]], CNVNATOR_STAT.out.pytor, [], [] )
    CNVNATOR_CALL ( [[:],[],[]], CNVNATOR_PARTITION.out.pytor, [], [] )
    CNVNATOR_CONVERT2VCF (CNVNATOR_CALL.out.tab)
}
