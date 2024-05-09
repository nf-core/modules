#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_INDEP } from '../../../../../modules/nf-core/plink/indep/main.nf'
include { PLINK_VCF } from '../../../../../modules/nf-core/plink/vcf/main.nf'

workflow test_plink_indep {

input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/genome/vcf/ped/justhusky_minimal.vcf.gz", checkIfExists: true)
    ]

    PLINK_VCF ( input )

    bed_ch = PLINK_VCF.out.bed
    bim_ch = PLINK_VCF.out.bim
    fam_ch = PLINK_VCF.out.fam

    ch_bed_bim_fam = bed_ch.join(bim_ch).join(fam_ch)

    ch_window_size = Channel.value(50)
    ch_variant_count = Channel.value(5)
    ch_variance_inflation_factor = Channel.value(1.5)

    PLINK_INDEP ( ch_bed_bim_fam, ch_window_size, ch_variant_count, ch_variance_inflation_factor )
}
