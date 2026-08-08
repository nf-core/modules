include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_NORMAL                                         } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TUMOR                                          } from '../../../modules/nf-core/samtools/index/main'
include { GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS as GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_NORMAL } from '../../../modules/nf-core/gridss/extractoverlappingfragments/main'
include { GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS as GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_TUMOR  } from '../../../modules/nf-core/gridss/extractoverlappingfragments/main'
include { GRIDSS_PREPROCESS as GRIDSS_PREPROCESS_NORMAL                                   } from '../../../modules/nf-core/gridss/preprocess/main'
include { GRIDSS_PREPROCESS as GRIDSS_PREPROCESS_TUMOR                                    } from '../../../modules/nf-core/gridss/preprocess/main'
include { GRIDSS_ASSEMBLE                                                                 } from '../../../modules/nf-core/gridss/assemble/main'
include { GRIDSS_CALL                                                                     } from '../../../modules/nf-core/gridss/call/main'
include { GRIDSS_SOMATICFILTER                                                            } from '../../../modules/nf-core/gridss/somaticfilter/main'


workflow BAM_TUMOR_NORMAL_SOMATIC_STRUCTURAL_VARIANT_CALLING_GRIDSS {

    take:
    ch_input_bam            // channel: [ val(meta), [ tumor_bam, normal_bam ], [tumor_bai, normal_bai] ]
    ch_fasta_fai_bwaindex   // channel: [ val(meta2), fasta, fai, bwa_index ]
    ch_pondir               // channel: [ val(meta3), [ pondir ] ]
    ch_gridss_config        // channel: [ val(meta4), [ gridss_config ] ]
    val_target_bed          // string:  target BED for panel/exome data, null to run untargeted

    main:

    def ch_input_normal = ch_input_bam.map { meta, bams, bais -> tuple(meta, bams[1], bais[1]) }
    def ch_input_tumor  = ch_input_bam.map { meta, bams, bais -> tuple(meta, bams[0], bais[0]) }

    // Targeted (panel/exome) data is subset to the fragments overlapping the target
    // regions before preprocessing, untargeted data goes straight to preprocessing
    if (val_target_bed) {

        def ch_target_bed = channel.value(tuple([ id:'target_bed' ], file(val_target_bed, checkIfExists: true)))

        GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_NORMAL ( ch_input_normal , ch_target_bed )
        GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_TUMOR  ( ch_input_tumor  , ch_target_bed )

        SAMTOOLS_INDEX_NORMAL ( GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_NORMAL.out.bam )
        SAMTOOLS_INDEX_TUMOR  ( GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_TUMOR.out.bam )

        ch_input_normal = GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_NORMAL.out.bam.join(SAMTOOLS_INDEX_NORMAL.out.index)
        ch_input_tumor  = GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS_TUMOR.out.bam.join(SAMTOOLS_INDEX_TUMOR.out.index)
    }

    GRIDSS_PREPROCESS_NORMAL ( ch_input_normal, ch_fasta_fai_bwaindex)
    GRIDSS_PREPROCESS_TUMOR  ( ch_input_tumor , ch_fasta_fai_bwaindex)

    def ch_input_assemble = ch_input_tumor
                            .join(ch_input_normal)
                            .join(GRIDSS_PREPROCESS_TUMOR.out.preprocess_dir)
                            .join(GRIDSS_PREPROCESS_NORMAL.out.preprocess_dir)
                            .map { meta, tbam, tbai, nbam, nbai, preprocess_tumor, preprocess_normal ->
                                tuple(meta, [tbam, nbam], [tbai, nbai], [preprocess_tumor, preprocess_normal])
                            }

    GRIDSS_ASSEMBLE ( ch_input_assemble, ch_fasta_fai_bwaindex, ch_gridss_config )

    def ch_input_call = ch_input_assemble.join(GRIDSS_ASSEMBLE.out.assemble_dir)

    GRIDSS_CALL ( ch_input_call, ch_fasta_fai_bwaindex, ch_gridss_config )

    GRIDSS_SOMATICFILTER ( GRIDSS_CALL.out.vcf, ch_pondir )


    emit:
    gridss_vcf                  = GRIDSS_CALL.out.vcf                      // channel: [ val(meta), join_call_vcf ]
    all_somatic_vcf             = GRIDSS_SOMATICFILTER.out.all_sv          // channel: [ val(meta), somatic_vcf ]
    high_confidence_somatic_vcf = GRIDSS_SOMATICFILTER.out.high_conf_sv    // channel: [ val(meta), high_conf_sv_vcf ]

}
