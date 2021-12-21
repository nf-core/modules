//
// Run SNPEFF to annotate VCF files
//

include { SNPEFF                                    } from '../../../modules/snpeff/main'
include { TABIX_BGZIPTABIX as ANNOTATION_BGZIPTABIX } from '../../../modules/tabix/bgziptabix/main'

workflow ANNOTATION_SNPEFF {
    take:
    vcf            // channel: [ val(meta), vcf ]
    snpeff_db      //   value: version of db to use
    snpeff_cache   //    path: path_to_snpeff_cache (optionnal)

    main:
    SNPEFF(vcf, snpeff_db, snpeff_cache)
    ANNOTATION_BGZIPTABIX(SNPEFF.out.vcf)
    ch_versions = SNPEFF.out.versions.first().mix(ANNOTATION_BGZIPTABIX.out.versions.first())

    emit:
    vcf_tbi  = ANNOTATION_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = SNPEFF.out.report                //    path: *.html
    versions = ch_versions                      //    path: versions.yml
}
