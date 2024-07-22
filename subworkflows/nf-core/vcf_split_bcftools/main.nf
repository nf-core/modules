include { BCFTOOLS_PLUGINSPLIT  } from '../../../modules/nf-core/bcftools/pluginsplit'

workflow VCF_SPLIT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr, tools], vcf ]

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_PLUGINSPLIT(ch_vcf, [], [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSPLIT.out.versions.first())

    ch_vcf_samples = BCFTOOLS_PLUGINSPLIT.out.vcf
        .transpose()
        .map{metaITC, vcf -> [metaITC + [id: vcf.getBaseName().tokenize(".")[0]], vcf]}

    ch_index_samples = BCFTOOLS_PLUGINSPLIT.out.merged_variants_index
        .transpose()
        .map{metaITC, index -> [metaITC + [id: index.getBaseName().tokenize(".")[0]], index]}

    ch_vcf_index_samples = ch_vcf_samples
        .join(ch_index_samples)

    emit:
    vcf_index      = ch_vcf_index_samples   // channel: [ [id, chr, tools], vcf, index ]
    versions       = ch_versions          // channel: [ versions.yml ]

}
