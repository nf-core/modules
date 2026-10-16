//
// Filter a structural variant VCF with SURVIVOR and bcftools
//

include { HTSLIB_BGZIPTABIX                            } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { HTSLIB_BGZIPTABIX as HTSLIB_BGZIPTABIX_BGZIP } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { SURVIVOR_FILTER                              } from '../../../modules/nf-core/survivor/filter/main'
include { BCFTOOLS_FILTER                              } from '../../../modules/nf-core/bcftools/filter/main'

// Please note that this subworkflow requires the options for SURVIVOR_FILTER defined in the nextflow.config of this directory
workflow VCF_FILTER_SURVIVOR_BCFTOOLS {

    take:
    ch_vcf        // channel: [ val(meta), path(vcf), path(bed) ] - pass [] for bed to filter genome-wide
    minsv         //   value: min SV size (-1 to disable)
    maxsv         //   value: max SV size (-1 to disable)
    minallelefreq //   value: min allele frequency (0-1)
    minnumreads   //   value: min number of reads supporting the SV (-1 to disable)

    main:
    // SURVIVOR can only read uncompressed VCFs
    HTSLIB_BGZIPTABIX(
        ch_vcf.map { meta, vcf, _bed -> [meta, vcf, [], []] },
        'decompress',
        false,
        '',
    )

    SURVIVOR_FILTER(
        HTSLIB_BGZIPTABIX.out.output
            .join(ch_vcf.map { meta, _vcf, bed -> [meta, bed] }, by: [0]),
        minsv,
        maxsv,
        minallelefreq,
        minnumreads,
    )

    // Filter expressions can be passed to bcftools filter with ext.args, e.g. '-e "QUAL<20"'
    BCFTOOLS_FILTER(
        SURVIVOR_FILTER.out.vcf.map { meta, vcf -> [meta, vcf, []] },
    )

    HTSLIB_BGZIPTABIX_BGZIP(
        BCFTOOLS_FILTER.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    ch_vcf_tbi = HTSLIB_BGZIPTABIX_BGZIP.out.output.join(HTSLIB_BGZIPTABIX_BGZIP.out.index, by: [0])

    emit:
    vcf     = HTSLIB_BGZIPTABIX_BGZIP.out.output // channel: [ val(meta), path(vcf.gz) ]
    tbi     = HTSLIB_BGZIPTABIX_BGZIP.out.index  // channel: [ val(meta), path(tbi) ]
    vcf_tbi = ch_vcf_tbi                         // channel: [ val(meta), path(vcf.gz), path(tbi) ]
}
