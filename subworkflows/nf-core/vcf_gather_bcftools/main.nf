include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_GATHER_BCFTOOLS {

    take:
    ch_vcfs             // channel: [ meta, vcf, tbi ]
    ch_scatter_output   // channel: [ meta, bed, gather_count ] => output from the scatter subworkflow, if you didn't use this subworkflow you can just use `[]` as bed since it isn't used
    val_common_meta     // string:  The name of the meta field that should become the new id
    val_sort            // boolean: Whether or not the output file should be sorted !! Add the config when using sort !!

    main:

    ch_versions = Channel.empty()

    ch_concat_input = ch_vcfs.join(ch_scatter_output)
        .map{ meta, vcf, tbi, bed, gather_count ->
            meta = val_common_meta ? meta + [id:meta[val_common_meta]] : meta
            [ groupKey(meta, gather_count), vcf, tbi ]
        }.groupTuple()

    BCFTOOLS_CONCAT ( ch_concat_input )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    if (val_sort) {
        BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

        ch_tabix_input = BCFTOOLS_SORT.out.vcf

    } else {
        ch_tabix_input = BCFTOOLS_CONCAT.out.vcf
    }

    TABIX_TABIX ( ch_tabix_input )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf      = ch_tabix_input        // channel: [ val(meta), [ vcf ] ]
    tbi      = TABIX_TABIX.out.tbi   // channel: [ val(meta), [ tbi ] ]

    versions = ch_versions           // channel: [ versions.yml ]
}
