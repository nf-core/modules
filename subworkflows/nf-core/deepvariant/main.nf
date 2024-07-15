
include { DEEPVARIANT_MAKE_EXAMPLES        } from '../../../modules/nf-core/deepvariant/makeexamples/main'
include { DEEPVARIANT_CALL_VARIANTS        } from '../../../modules/nf-core/deepvariant/callvariants/main'
include { DEEPVARIANT_POSTPROCESS_VARIANTS } from '../../../modules/nf-core/deepvariant/postprocessvariants/main'

workflow DEEPVARIANT {
    take:
    ch_input            // channel: [ val(meta), path(input), path(index), path(intervals)]
    ch_fasta            // channel: [ val(meta2), path(fasta) ]
    ch_fai              // channel: [ val(meta3), path(fail) ]
    ch_gzi              // channel: [ val(meta4), path(gzi) ]
    ch_model_type       // channel: val("wgs" | "wes")

    main:

    ch_versions = Channel.empty()

    DEEPVARIANT_MAKE_EXAMPLES(ch_input, ch_fasta, ch_fai, ch_gzi)
    ch_versions = ch_versions.mix(DEEPVARIANT_MAKE_EXAMPLES.out.versions.first())

    DEEPVARIANT_CALL_VARIANTS(DEEPVARIANT_MAKE_EXAMPLES.out.examples, ch_model_type)
    ch_versions = ch_versions.mix(DEEPVARIANT_CALL_VARIANTS.out.versions.first())
    
    // Input to postprocessing step: each item. for each unique (meta, intervals),
    // should contain one set of gvcf records and one variant call tfrecord.
    ch_postproc_input = DEEPVARIANT_MAKE_EXAMPLES.out.gvcf.join(
        DEEPVARIANT_CALL_VARIANTS.out.call_variants_tfrecords,
        by: [0, 1],
        failOnMismatch: true
    )
    DEEPVARIANT_POSTPROCESS_VARIANTS(
        ch_postproc_input,
        ch_fasta,
        ch_fai,
        ch_gzi
    )    
    ch_versions = ch_versions.mix(DEEPVARIANT_POSTPROCESS_VARIANTS.out.versions.first())


    emit:
    vcf         = DEEPVARIANT_POSTPROCESS_VARIANTS.out.vcf
    vcf_tbi     = DEEPVARIANT_POSTPROCESS_VARIANTS.out.vcf_tbi
    gvcf        = DEEPVARIANT_POSTPROCESS_VARIANTS.out.gvcf
    gvcf_tbi    = DEEPVARIANT_POSTPROCESS_VARIANTS.out.gvcf_tbi
    versions    = ch_versions
}
