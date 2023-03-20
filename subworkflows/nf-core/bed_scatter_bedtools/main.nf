include { BEDTOOLS_SPLIT     } from '../../../modules/nf-core/bedtools/split/main'

workflow BED_SCATTER_BEDTOOLS {

    // IMPORTANT: Add the configuration in nextflow.config file to your modules.config

    take:
    ch_bed // channel: [ val(meta), path(bed), val(scatter_count) ]

    main:

    ch_versions = Channel.empty()

    ch_bedtools_split = ch_bed
        .map(
            { meta, bed, scatter_count ->
                meta = meta + [subwf_scatter_count:scatter_count]
                [ meta, bed ]
            }
        )

    BEDTOOLS_SPLIT(
        ch_bedtools_split
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions.first())

    ch_scattered_beds = BEDTOOLS_SPLIT.out.beds
        .map(
            { meta, beds ->
                // Checks if the scatter count corresponds to the amount of files created. (This doesn't match in some edge cases)
                scatter_count = beds instanceof Path ? 1 : beds.size()
                meta.remove("subwf_scatter_count")
                [ meta, beds, scatter_count ]
            }
        )
        .transpose()

    emit:
    scattered_beds  = ch_scattered_beds // channel: [ val(meta), path(bed), val(scatter_count) ]

    versions        = ch_versions       // channel: [ path(versions.yml) ]
}

