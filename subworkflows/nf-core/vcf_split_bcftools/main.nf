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

    ch_tbi_samples = BCFTOOLS_PLUGINSPLIT.out.tbi
        .transpose()
        .map{metaITC, tbi -> [metaITC + [id: tbi.getBaseName().tokenize(".")[0]], tbi]}

    ch_vcf_tbi_samples = ch_vcf_samples
        .join(ch_tbi_samples)

    emit:
    vcf_tbi        = ch_vcf_tbi_samples   // channel: [ [id, chr, tools], vcf, index ]
    versions       = ch_versions          // channel: [ versions.yml ]

}
