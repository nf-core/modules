#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLINK_EXCLUDE } from '../../../../../modules/nf-core/plink/exclude/main.nf'
include { PLINK_VCF     } from '../../../../../modules/nf-core/plink/vcf/main.nf'
include { GAWK          } from '../../../../../modules/nf-core/gawk/main.nf'

workflow test_plink_exclude {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['syntheticvcf_short_vcf_gz'], checkIfExists: true)
    ]

    PLINK_VCF ( input )
    bed_ch = PLINK_VCF.out.bed
    bim_ch = PLINK_VCF.out.bim
    fam_ch = PLINK_VCF.out.fam

    ch_bed_bim_fam = bed_ch.join(bim_ch).join(fam_ch)

    ch_program_file = Channel.value('BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA")  {print $2}').collectFile(name:"program.txt")
    GAWK( bim_ch , ch_program_file)
    ch_variants = GAWK.out.output

    ch_bed_bim_fam_variants = ch_bed_bim_fam.join(ch_variants)

    PLINK_EXCLUDE ( ch_bed_bim_fam_variants )
}
