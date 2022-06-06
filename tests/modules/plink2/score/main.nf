#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK2_VCF } from '../../../../modules/plink2/vcf/main.nf'
include { PLINK2_SCORE } from '../../../../modules/plink2/score/main.nf'

workflow test_plink2_score {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true)
    ]
    PLINK2_VCF ( input )

    scorefile = file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_score'], checkIfExists: true)

    PLINK2_VCF.out.pgen
        .concat(PLINK2_VCF.out.psam, PLINK2_VCF.out.pvar)
        .groupTuple()
        .map { it.flatten() }
        .set { ch_target_genome }

    PLINK2_SCORE ( ch_target_genome, scorefile )
}
