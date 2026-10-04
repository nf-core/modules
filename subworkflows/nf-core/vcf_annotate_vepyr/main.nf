//
// Optionally download a cache and normalise a VCF, then annotate it with vepyr
//

include { BCFTOOLS_NORM  } from '../../../modules/nf-core/bcftools/norm/main'
include { HUGGINGFACE_DOWNLOAD } from '../../../modules/nf-core/huggingface/download/main'
include { VEPYR_ANNOTATE } from '../../../modules/nf-core/vepyr/annotate/main'

workflow VCF_ANNOTATE_VEPYR {
    take:
    ch_vcf            // channel: [ val(meta), path(vcf), path(tbi) ] (tbi optional, pass [])
    ch_cache          // channel: [ val(meta2), path(cache) ] (value channel; pass [[], []] when val_hf_repo is set)
    ch_fasta          // channel: [ val(meta3), path(fasta), path(fai), path(gzi) ] (value channel; gzi only for a bgzip FASTA, else [])
    val_cache_version // value:   integer, or [] to skip the cache version check
    ch_plugin_cache   // channel: [ val(meta4), path(plugin_cache) ] (optional, pass [ [], [] ])
    val_normalize     // value:   boolean, run BCFTOOLS_NORM before annotation
    val_hf_repo       // value:   Hugging Face dataset ID, or '' / [] to use ch_cache

    main:
    def ch_cache_value = ch_cache instanceof List ? channel.value(ch_cache) : ch_cache

    if (val_hf_repo) {
        def hf_repo = val_hf_repo.startsWith('datasets/') ? val_hf_repo : "datasets/${val_hf_repo}"
        // A value channel shares one cache download across all input samples.
        HUGGINGFACE_DOWNLOAD(channel.value([
            [id: val_hf_repo.tokenize('/').last()],
            hf_repo,
            [],
        ]))
        ch_annotate_cache = HUGGINGFACE_DOWNLOAD.out.output
            .combine(ch_cache_value.map { _meta, cache ->
                if (cache) {
                    error("VCF_ANNOTATE_VEPYR: set either ch_cache or val_hf_repo, not both.")
                }
                true
            })
            .map { meta, cache, _no_local_cache -> [meta, cache] }
            .first()
    }
    else {
        ch_annotate_cache = ch_cache_value.map { meta, cache ->
            if (!cache) {
                error("VCF_ANNOTATE_VEPYR: no cache given, set ch_cache or val_hf_repo.")
            }
            [meta, cache]
        }
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
