// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { NGSCHECKMATE_NCM } from '../../../modules/nf-core/ngscheckmate/ncm/main'
// please note this subworkflow requires the options for bcltools_mpileup that are included in the nextflow.config
workflow BAM_NGSCHECKMATE {
    take:
    ch_input // channel: [ val(meta1), bam/cram ]
    ch_snp_bed // channel: [ val(meta2), bed ]
    ch_fasta // channel: [ val(meta3), fasta, fai ]

    main:
    ch_input_bed = ch_input
        .combine(ch_snp_bed)
        .map { input_meta, input_file, _bed_meta, bed_file ->
            [input_meta, input_file, bed_file, []]
        }

    BCFTOOLS_MPILEUP(ch_input_bed, ch_fasta.collect(), false)

    BCFTOOLS_MPILEUP.out.vcf.map { _meta, vcf -> vcf }.collect().map { files -> [files] }.set { ch_collected_vcfs }

    ch_snp_bed
        .map { meta, _bed -> meta }
        .combine(ch_collected_vcfs)
        .set { ch_vcfs }

    NGSCHECKMATE_NCM(ch_vcfs, ch_snp_bed, ch_fasta)

    emit:
    corr_matrix = NGSCHECKMATE_NCM.out.corr_matrix // channel: [ meta, corr_matrix ]
    matched     = NGSCHECKMATE_NCM.out.matched // channel: [ meta, matched ]
    all         = NGSCHECKMATE_NCM.out.all // channel: [ meta, all ]
    vcf         = BCFTOOLS_MPILEUP.out.vcf // channel: [ meta, vcf ]
    pdf         = NGSCHECKMATE_NCM.out.pdf // channel: [ meta, pdf ]
}
