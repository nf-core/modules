//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    ch_bedgraph // channel: [ val(meta), path(bedgraph) ]
    ch_sizes    // channel: [ path (chrom.sizes)]

    main:

    ch_versions = Channel.empty()

    //
    // Clip bedGraph file
    //
    UCSC_BEDCLIP ( ch_bedgraph, ch_sizes )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())

    //
    // Convert bedGraph to bigWig
    //
    UCSC_BEDGRAPHTOBIGWIG ( UCSC_BEDCLIP.out.ch_bedgraph, ch_sizes )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), path(bigwig) ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), path(bedgraph) ]

    versions = ch_versions                      // channel: [ path(versions.yml) ]
}
