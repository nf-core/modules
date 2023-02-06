#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ADMIXTURE } from '../../../../modules/nf-core/admixture/main.nf'
include { PLINK_VCF } from '../../../../modules/nf-core/plink/vcf/main.nf'


workflow test_admixture {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/illumina/vcf/test_annotate.vcf.gz", checkIfExists: true)
    ]
    PLINK_VCF ( input )

    bed_ch = PLINK_VCF.out.bed
    bim_ch = PLINK_VCF.out.bim
    fam_ch = PLINK_VCF.out.fam

    ch_bed_bim_fam = bed_ch.join(bim_ch).join(fam_ch)

    ch_value_K = Channel.of("3")

    ADMIXTURE ( ch_bed_bim_fam, ch_value_K )

}
