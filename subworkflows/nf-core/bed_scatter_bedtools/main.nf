include { BEDTOOLS_SPLIT } from '../../../modules/nf-core/bedtools/split'

workflow BED_SCATTER_BEDTOOLS {
    take:
    ch_bed // channel: [ val(meta), path(bed), val(scatter_count) ]

    main:
    BEDTOOLS_SPLIT(ch_bed)

    // Checks if the scatter count corresponds to the amount of files created. (This doesn't match in some edge cases)
    ch_scattered_beds = BEDTOOLS_SPLIT.out.beds
        .map { meta, beds ->
            def scatter_count = beds instanceof Path ? 1 : beds.size()
            [meta, beds, scatter_count]
        }
        .transpose()

    emit:
    scattered_beds = ch_scattered_beds // channel: [ val(meta), path(bed), val(scatter_count) ]
}
