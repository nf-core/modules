
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
    // Assign a unique ID to each input element, so the input can be processed
    // as a single unique item by all subworkflow stages. This is necessary to merge
    // the GVCF data and called variants in the postprocessing step.
    ch_input_with_uid = ch_input.map {
        meta, input, index, intervals -> [ ['deepvariant_id': UUID.randomUUID()] + meta, input, index, intervals]
    }

    DEEPVARIANT_MAKEEXAMPLES(ch_input_with_uid, ch_fasta, ch_fai, ch_gzi)
    ch_versions = ch_versions.mix(DEEPVARIANT_MAKEEXAMPLES.out.versions.first())

    DEEPVARIANT_CALLVARIANTS(DEEPVARIANT_MAKEEXAMPLES.out.examples, ch_model_type)
    ch_versions = ch_versions.mix(DEEPVARIANT_CALLVARIANTS.out.versions.first())
    
    // Input to postprocessing step needs both the gvcfs from MAKEEXAMPLES and the variant
    // calls from CALLVARIANTS. Uniqueness of the joining key is guaranteed by the deepvariant_id
    // element added to meta.
    ch_postproc_input = DEEPVARIANT_CALLVARIANTS.out.call_variants_tfrecords.join(
        DEEPVARIANT_MAKEEXAMPLES.out.gvcf,
        failOnMismatch: true
    )
    DEEPVARIANT_POSTPROCESSVARIANTS(
        ch_postproc_input,
        ch_fasta,
        ch_fai,
        ch_gzi
    )

    // Create output channels without deepvariant_id
    def remove_deepvariant_id = {
        newMeta = it[0];
        newMeta.remove("deepvariant_id");
        [newMeta] + it.drop(1)
    }
    ch_versions = ch_versions.mix(DEEPVARIANT_POSTPROCESSVARIANTS.out.versions.first())

    emit:
    vcf         = DEEPVARIANT_POSTPROCESSVARIANTS.out.vcf.map(remove_deepvariant_id)
    vcf_tbi     = DEEPVARIANT_POSTPROCESSVARIANTS.out.vcf_tbi.map(remove_deepvariant_id)
    gvcf        = DEEPVARIANT_POSTPROCESSVARIANTS.out.gvcf.map(remove_deepvariant_id)
    gvcf_tbi    = DEEPVARIANT_POSTPROCESSVARIANTS.out.gvcf_tbi.map(remove_deepvariant_id)
    versions    = ch_versions
}
