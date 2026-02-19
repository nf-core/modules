include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT   } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX     } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_GATHER_BCFTOOLS {
    take:
    ch_vcfs           // channel: [ meta, vcf, index, count ]
    arr_common_meta   // array: The name of the meta fields that should be used for grouping
    val_sort          // boolean: Whether or not the output file should be sorted !! Add the config when using sort !!

    main:

    // Check if arr_common_meta is an array or list
    if (!(arr_common_meta instanceof List || arr_common_meta instanceof Collection)) {
        error("ERROR: arr_common_meta should be an array or list, got ${arr_common_meta.getClass()}")
    }

    ch_concat_input = ch_vcfs
        .map { meta, vcf, index, count ->
            def missingKeys = arr_common_meta.findAll { key -> !(key in meta) }
            if (missingKeys) {
                error("ERROR: Keys ${missingKeys} from arr_common_meta not found in meta. Available keys: ${meta.keySet()}")
            }
            def newMeta = arr_common_meta ?
                arr_common_meta.collectEntries { key -> [(key): meta[key]] } :
                meta
            [groupKey(newMeta, count), meta, vcf, index]
        }
        .groupTuple()
        .ifEmpty {
            error("ERROR: grouping operation resulted in an empty channel.")
        }
        .branch { key, meta, vcf, index ->
            def cleanedMetas = meta.collect { m ->
                m.findAll { k, _v -> !(k in arr_common_meta) }
            }
            def newMeta = arr_common_meta ? key.target + [metas: cleanedMetas] : meta[0]
            def out_tuple = [newMeta, vcf, index]

            one: vcf.size() == 1
                return out_tuple
            more: vcf.size() > 1
                return out_tuple
        }

    // Concatenate vcf with more than one record
    BCFTOOLS_CONCAT(ch_concat_input.more)

    ch_vcf_concat = ch_concat_input.one
        .map{ meta, vcf, _index -> [meta, vcf.get(0)] }
        .mix(BCFTOOLS_CONCAT.out.vcf)

    if (val_sort) {
        BCFTOOLS_SORT(ch_vcf_concat)
        ch_tabix_input = BCFTOOLS_SORT.out.vcf
    } else {
        ch_tabix_input = ch_vcf_concat
    }

    TABIX_TABIX(ch_tabix_input)

    ch_vcf_index = ch_tabix_input
        .join(TABIX_TABIX.out.index)

    emit:
    vcf       = ch_vcf_concat         // channel: [ val(meta), [ vcf ] ]
    index     = TABIX_TABIX.out.index // channel: [ val(meta), [ tbi or csi ] ]
    vcf_index = ch_vcf_index          // channel: [ val(meta), [ vcf ], [ tbi or csi ] ]
}
