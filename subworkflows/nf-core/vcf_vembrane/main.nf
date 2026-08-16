include { TABIX_BGZIPTABIX                       } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { VEMBRANE_FILTER                        } from '../../../modules/nf-core/vembrane/filter/main'
include { VEMBRANE_SORT as VEMBRANE_SORT_DEFAULT } from '../../../modules/nf-core/vembrane/sort/main'
include { VEMBRANE_SORT as VEMBRANE_SORT_VCF     } from '../../../modules/nf-core/vembrane/sort/main'
include { VEMBRANE_TABLE                         } from '../../../modules/nf-core/vembrane/table/main'

workflow VCF_VEMBRANE {
    take:
    ch_vcf            // channel: [ val(meta), path(vcf) ]
    filter_expression // string:  [mandatory] expression for vembrane filter
    sort_expression   // string:  [mandatory] expression for vembrane sort
    table_expression  // string:  [mandatory] expression for vembrane table

    main:
    ch_versions = Channel.empty()

    // Start with original VCF channel
    ch_vcf_current = ch_vcf

    // Conditionally FILTER VCF by vembrane
    if (filter_expression && filter_expression.toString().trim() != "" && filter_expression.toString().trim() != "null") {
        VEMBRANE_FILTER(ch_vcf_current, filter_expression)
        ch_vcf_current = VEMBRANE_FILTER.out.filtered_variant
        ch_versions = ch_versions.mix(VEMBRANE_FILTER.out.versions)
    }

    // Conditionally SORT VCF by vembrane
    if (sort_expression && sort_expression.toString().trim() != "" && sort_expression.toString().trim() != "null") {
        VEMBRANE_SORT_VCF(ch_vcf_current, sort_expression)
        ch_vcf_current = VEMBRANE_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(VEMBRANE_SORT_VCF.out.versions)
    }

    // Always sort by 'CHROM, POS' before index and output compression
    VEMBRANE_SORT_DEFAULT(ch_vcf_current, "CHROM, POS")
    ch_vcf_current = VEMBRANE_SORT_DEFAULT.out.vcf
    ch_versions = ch_versions.mix(VEMBRANE_SORT_DEFAULT.out.versions)

    // Index VCF.GZ (always performed on final VCF)
    TABIX_BGZIPTABIX(ch_vcf_current)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    // Conditionally create TSV table
    if (table_expression && table_expression.toString().trim() != "" && table_expression.toString().trim() != "null") {
        VEMBRANE_TABLE(ch_vcf_current, table_expression)
        tsv_ch = VEMBRANE_TABLE.out.table
        ch_versions = ch_versions.mix(VEMBRANE_TABLE.out.versions)
    }
    else {
        // Create empty channel when table is not generated
        tsv_ch = Channel.empty()
    }
    ch_vcf = TABIX_BGZIPTABIX.out.gz_tbi

    emit:
    vcf      = ch_vcf_current // channel: [ val(meta), path(vcf_file), path(tbi_file) ]
    tsv      = tsv_ch         // channel: [ val(meta), path(table_file) ] or empty

    versions = ch_versions    // channel: [ path(versions.yml) ]
}
