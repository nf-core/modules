include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_GATHER_BCFTOOLS {

    take:
    ch_vcfs             // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_scatter_count    // channel: [ val(meta), val(gather_count) ] => The scatter count per group of input files
    ch_bed              // channel: [ path(bed) ] => The BED file to be used by bcftools concat
    val_common_meta     // string:  The name of the meta field that should become the new id
    val_sort            // boolean: Whether or not the output file should be sorted !! Add the config when using sort !!

    main:

    ch_versions = Channel.empty()

    ch_concat_input = ch_vcfs.join(ch_scatter_count)
        .map{ meta, vcf, tbi, gather_count ->
            new_meta = val_common_meta ? meta + [id:meta[val_common_meta]] : meta
            [ groupKey(new_meta, gather_count), vcf, tbi ]
        }.groupTuple()

    BCFTOOLS_CONCAT ( ch_concat_input, ch_bed )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    if (val_sort) {
        BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

        ch_tabix_input = BCFTOOLS_SORT.out.vcf

    } else {
        ch_tabix_input = BCFTOOLS_CONCAT.out.vcf
    }

    TABIX_TABIX ( ch_tabix_input )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    emit:
    vcf      = ch_tabix_input        // channel: [ val(meta), path(vcf) ]
    tbi      = TABIX_TABIX.out.tbi   // channel: [ val(meta), path(tbi) ]

    versions = ch_versions           // channel: [ versions.yml ]
}

