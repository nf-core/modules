//
// Optionally download a cache and normalise a VCF, then annotate it with vepyr
//

include { BCFTOOLS_NORM  } from '../../../modules/nf-core/bcftools/norm/main'
include { HUGGINGFACE_DOWNLOAD } from '../../../modules/nf-core/huggingface/download/main'
include { VEPYR_ANNOTATE } from '../../../modules/nf-core/vepyr/annotate/main'

workflow VCF_ANNOTATE_VEPYR {
    take:
    ch_vcf            // channel: [ val(meta), path(vcf), path(tbi) ] (tbi optional, pass [])
    ch_cache          // channel: [ val(meta2), path(cache) ] (unused when val_hf_repo is set, pass [[], []])
    ch_fasta          // channel: [ val(meta3), path(fasta), path(fai), path(gzi) ] (gzi only for a bgzip FASTA, else [])
    val_cache_version // value:   integer, or [] to skip the cache version check
    ch_plugin_cache   // channel: [ val(meta4), path(plugin_cache) ] (optional, pass [ [], [] ])
    val_normalize     // value:   boolean, run BCFTOOLS_NORM before annotation (recommended: true)
    val_hf_repo       // value:   Hugging Face cache repository ID, or '' / [] to use ch_cache

    main:
    ch_annotate_cache = ch_cache
    if (val_hf_repo) {
        // A value channel shares one cache download across all input samples.
        HUGGINGFACE_DOWNLOAD(channel.value([
            [id: val_hf_repo.tokenize('/').last()],
            val_hf_repo,
            [],
        ]))
        ch_annotate_cache = HUGGINGFACE_DOWNLOAD.out.output
    }

    ch_annotate_input = ch_vcf

    if (val_normalize) {
        BCFTOOLS_NORM(
            ch_vcf,
            ch_fasta.map { meta, fasta, _fai, _gzi -> [meta, fasta] },
        )

        // Without --write-index bcftools/norm emits no index; VEPYR_ANNOTATE then
        // builds one.
        ch_annotate_input = BCFTOOLS_NORM.out.vcf
            .join(BCFTOOLS_NORM.out.index, failOnDuplicate: true, remainder: true)
            .map { meta, vcf, index -> [meta, vcf, index ?: []] }
    }

    VEPYR_ANNOTATE(
        ch_annotate_input,
        ch_annotate_cache,
        ch_fasta,
        val_cache_version,
        ch_plugin_cache,
    )

    ch_vcf_tbi = VEPYR_ANNOTATE.out.vcf.join(VEPYR_ANNOTATE.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    emit:
    vcf_tbi = ch_vcf_tbi // channel: [ val(meta), path(vcf), path(tbi) ]
}
