include { ENSEMBLVEP_FILTERVEP } from '../../../modules/nf-core/ensemblvep/filtervep/main'
include { BCFTOOLS_VIEW        } from '../../../modules/nf-core/bcftools/view/main'
include { TABIX_BGZIPTABIX     } from '../../../modules/nf-core/tabix/bgziptabix/main'
// Please note this subworkflow requires the options for bcftools_view that are included in the nextflow.config
workflow VCF_FILTER_BCFTOOLS_ENSEMBLVEP {

    take:
    ch_vcf                     // channel: [ val(meta), path(vcf) ]
    ch_filter_vep_feature_file // channel: [ val(meta), path(txt) ]
    filter_with_bcftools       //    bool: should bcftools view be run
    filter_with_filter_vep     //    bool: should filter_vep be run

    main:
    ch_versions = Channel.empty()
    ch_tbi = Channel.empty()

    // Since bcftools is likely much faster than filter_vep,
    // we run it first to reduce the number of variants that filter_vep has to process.
    if(filter_with_bcftools) {

        BCFTOOLS_VIEW (
            ch_vcf.map { meta, vcf -> [ meta, vcf, [] ] },
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        ch_vcf = BCFTOOLS_VIEW.out.vcf
        ch_tbi = BCFTOOLS_VIEW.out.tbi
    }

    if(filter_with_filter_vep) {

        ENSEMBLVEP_FILTERVEP (
            ch_vcf,
            ch_filter_vep_feature_file.map { _meta, file -> file }
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_FILTERVEP.out.versions)

        TABIX_BGZIPTABIX (
            ENSEMBLVEP_FILTERVEP.out.output
        )
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

        ch_vcf = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, _tbi -> [ meta, vcf ] }
        ch_tbi = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _vcf, tbi -> [ meta, tbi ] }
    }

    emit:
    vcf      = ch_vcf      // channel: [ val(meta), path(vcf) ]
    tbi      = ch_tbi      // channel: [ val(meta), path(tbi) ]
    versions = ch_versions // channel: [ path(versions.yml) ]
}
