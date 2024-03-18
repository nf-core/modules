#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_RECODE } from '../../../../../modules/nf-core/plink/recode/main.nf'
include { PLINK_VCF } from '../../../../../modules/nf-core/plink/vcf/main.nf'

workflow test_plink_recode {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true)
    ]
    PLINK_VCF ( input )

    bed_ch = PLINK_VCF.out.bed
    bim_ch = PLINK_VCF.out.bim
    fam_ch = PLINK_VCF.out.fam

    ch_bed_bim_fam = bed_ch.join(bim_ch).join(fam_ch)

    PLINK_RECODE ( ch_bed_bim_fam )
}
