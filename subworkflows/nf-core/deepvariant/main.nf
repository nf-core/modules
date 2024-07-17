
include { DEEPVARIANT_MAKEEXAMPLES        } from '../../../modules/nf-core/deepvariant/makeexamples/main'
include { DEEPVARIANT_CALLVARIANTS        } from '../../../modules/nf-core/deepvariant/callvariants/main'
include { DEEPVARIANT_POSTPROCESSVARIANTS } from '../../../modules/nf-core/deepvariant/postprocessvariants/main'

workflow DEEPVARIANT {
    take:
    ch_input            // channel: [ val(meta), path(input), path(index), path(intervals)]
    ch_fasta            // channel: [ val(meta2), path(fasta) ]
    ch_fai              // channel: [ val(meta3), path(fail) ]
    ch_gzi              // channel: [ val(meta4), path(gzi) ]
    ch_model_type       // channel: val("wgs" | "wes")

    main:

    ch_versions = Channel.empty()

    DEEPVARIANT_MAKEEXAMPLES(ch_input, ch_fasta, ch_fai, ch_gzi)
    ch_versions = ch_versions.mix(DEEPVARIANT_MAKEEXAMPLES.out.versions.first())

    DEEPVARIANT_CALLVARIANTS(DEEPVARIANT_MAKEEXAMPLES.out.examples, ch_model_type)
    ch_versions = ch_versions.mix(DEEPVARIANT_CALLVARIANTS.out.versions.first())
    
    // Input to postprocessing step needs both the gvcfs from MAKEEXAMPLES and the variant
    // calls from CALLVARIANTS. Joining on the tuple element make_examples_id.
    ch_postproc_input = DEEPVARIANT_CALLVARIANTS.out.call_variants_tfrecords.join(
        DEEPVARIANT_MAKEEXAMPLES.out.gvcf,
        by: [0, 1], // Join on [0] is sufficient, but we want to collapse meta as well.
        failOnMismatch: true
    ).map { it.drop(1) } // Drop the unique ID, which shouldn't be input into POSTPROCESSVARIANTS
    DEEPVARIANT_POSTPROCESSVARIANTS(
        ch_postproc_input,
        ch_fasta,
        ch_fai,
        ch_gzi
    )

    ch_versions = ch_versions.mix(DEEPVARIANT_POSTPROCESSVARIANTS.out.versions.first())

    emit:
    vcf         = DEEPVARIANT_POSTPROCESSVARIANTS.out.vcf
    vcf_tbi     = DEEPVARIANT_POSTPROCESSVARIANTS.out.vcf_tbi
    gvcf        = DEEPVARIANT_POSTPROCESSVARIANTS.out.gvcf
    gvcf_tbi    = DEEPVARIANT_POSTPROCESSVARIANTS.out.gvcf_tbi
    versions    = ch_versions
}
