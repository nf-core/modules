include { BEDTOOLS_SPLIT     } from '../../../modules/nf-core/bedtools/split/main'

workflow BED_SCATTER_BEDTOOLS {

    take:
    ch_bed // channel: [ val(meta), path(bed), val(scatter_count) ]

    main:

    ch_versions = Channel.empty()

    BEDTOOLS_SPLIT(
        ch_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    ch_scattered_beds = BEDTOOLS_SPLIT.out.beds
        .map(
            { meta, beds ->
                // Checks if the scatter count corresponds to the amount of files created. (This doesn't match in some edge cases)
                scatter_count = beds instanceof Path ? 1 : beds.size()
                [ meta, beds, scatter_count ]
            }
        )
        .transpose()

    emit:
    scattered_beds  = ch_scattered_beds // channel: [ val(meta), path(bed), val(scatter_count) ]

    versions        = ch_versions       // channel: [ path(versions.yml) ]
}

