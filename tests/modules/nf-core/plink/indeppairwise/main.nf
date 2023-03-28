#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_INDEPPAIRWISE } from '../../../../../modules/nf-core/plink/indeppairwise/main.nf'
include { PLINK_VCF           } from '../../../../../modules/nf-core/plink/vcf/main.nf'

workflow test_plink_indeppairwise {

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
    ch_r2_threshold = Channel.value(0.2)

    PLINK_INDEPPAIRWISE ( ch_bed_bim_fam, ch_window_size, ch_variant_count, ch_r2_threshold )
}
