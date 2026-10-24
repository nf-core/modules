// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

//
// Filter a structural variant VCF with SURVIVOR and/or bcftools
//

include { BCFTOOLS_FILTER                                   } from '../../../modules/nf-core/bcftools/filter/main'
include { HTSLIB_BGZIPTABIX as HTSLIB_BGZIPTABIX_DECOMPRESS } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { HTSLIB_BGZIPTABIX as HTSLIB_BGZIPTABIX_COMPRESS   } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { SURVIVOR_FILTER                                   } from '../../../modules/nf-core/survivor/filter/main'

// Please note that this subworkflow requires the options for SURVIVOR_FILTER defined in the nextflow.config of this directory
workflow VCF_FILTER_SURVIVOR_BCFTOOLS {

    take:
    ch_vcf           // channel: [ val(meta), path(vcf), path(bed) ] - pass [] for bed to filter genome-wide
    minsv            //   value: min SV size (-1 to disable)
    maxsv            //   value: max SV size (-1 to disable)
    minallelefreq    //   value: min allele frequency (0-1)
    minnumreads      //   value: min number of reads supporting the SV (-1 to disable)
    val_tools_to_use //   value: a list of tools to use, options are: ["survivor", "bcftools"]

    main:
    def allowed_tools = ['survivor', 'bcftools']
    if (!(val_tools_to_use instanceof Collection) || val_tools_to_use.isEmpty() || !val_tools_to_use.every { tool -> tool in allowed_tools }) {
        error("VCF_FILTER_SURVIVOR_BCFTOOLS: `val_tools_to_use` must be a list of one or more of ${allowed_tools.join(', ')}, got: ${val_tools_to_use}")
    }

    def ch_filtered_vcf = ch_vcf.map { meta, vcf, _bed -> [meta, vcf] }

    // SURVIVOR can only read uncompressed VCFs
    if ('survivor' in val_tools_to_use) {
        HTSLIB_BGZIPTABIX_DECOMPRESS(
            ch_vcf.map { meta, vcf, _bed -> [meta, vcf, [], []] },
            'decompress',
            false,
            'vcf',
        )

        SURVIVOR_FILTER(
            HTSLIB_BGZIPTABIX_DECOMPRESS.out.output
                .join(ch_vcf.map { meta, _vcf, bed -> [meta, bed] }, failOnDuplicate: true, failOnMismatch: true),
            minsv,
            maxsv,
            minallelefreq,
            minnumreads,
        )

        ch_filtered_vcf = SURVIVOR_FILTER.out.vcf
    }

    // Filter expressions can be passed to bcftools filter with ext.args, e.g. '-e "QUAL<20"'
    if ('bcftools' in val_tools_to_use) {
        BCFTOOLS_FILTER(
            ch_filtered_vcf.map { meta, vcf -> [meta, vcf, []] },
        )

        ch_filtered_vcf = BCFTOOLS_FILTER.out.vcf
    }

    HTSLIB_BGZIPTABIX_COMPRESS(
        ch_filtered_vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    ch_vcf_tbi = HTSLIB_BGZIPTABIX_COMPRESS.out.output.join(HTSLIB_BGZIPTABIX_COMPRESS.out.index, failOnDuplicate: true, failOnMismatch: true)

    emit:
    vcf     = HTSLIB_BGZIPTABIX_COMPRESS.out.output // channel: [ val(meta), path(vcf.gz) ]
    tbi     = HTSLIB_BGZIPTABIX_COMPRESS.out.index  // channel: [ val(meta), path(tbi) ]
    vcf_tbi = ch_vcf_tbi                            // channel: [ val(meta), path(vcf.gz), path(tbi) ]
}
