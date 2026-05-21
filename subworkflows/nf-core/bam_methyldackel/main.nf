include { METHYLDACKEL_EXTRACT                          } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS                            } from '../../../modules/nf-core/methyldackel/mbias/main'

workflow BAM_METHYLDACKEL {

    take:
    ch_alignment            // channel: [ val(meta), [ bam ] ] ## BAM from alignment
    ch_alignment_index      // channel: [ val(meta), [ bai ] ] ## BAI from alignment
    ch_fasta                // channel: [ val(meta), [ fasta ] ]
    ch_fasta_index          // channel: [ val(meta), [ fasta index ] ]

    main:
    ch_methydackel_extract_bedgraph  = channel.empty()
    ch_methydackel_extract_methylkit = channel.empty()
    ch_methydackel_mbias             = channel.empty()
    ch_multiqc_files                 = channel.empty()
    ch_versions                      = channel.empty()

    /*
     * Extract per-base methylation and plot methylation bias
     */

    METHYLDACKEL_EXTRACT (
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_fasta_index
    )
    ch_methydackel_extract_bedgraph  = METHYLDACKEL_EXTRACT.out.bedgraph
    ch_methydackel_extract_methylkit = METHYLDACKEL_EXTRACT.out.methylkit
    ch_versions                      = ch_versions.mix(METHYLDACKEL_EXTRACT.out.versions)

    METHYLDACKEL_MBIAS (
        ch_alignment,
        ch_alignment_index,
        ch_fasta,
        ch_fasta_index
    )
    ch_methydackel_mbias = METHYLDACKEL_MBIAS.out.txt
    ch_versions          = ch_versions.mix(METHYLDACKEL_MBIAS.out.versions)

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_methydackel_extract_bedgraph.collect{ _meta, bedgraph -> bedgraph  }
                        .mix(ch_methydackel_extract_methylkit.collect{ _meta, methylkit -> methylkit })
                        .mix(ch_methydackel_mbias.collect{ _meta, txt -> txt  })

    emit:
    methydackel_extract_bedgraph  = ch_methydackel_extract_bedgraph  // channel: [ val(meta), [ bedgraph ]  ]
    methydackel_extract_methylkit = ch_methydackel_extract_methylkit // channel: [ val(meta), [ methylkit ] ]
    methydackel_mbias             = ch_methydackel_mbias             // channel: [ val(meta), [ mbias ]     ]
    multiqc                       = ch_multiqc_files                 // channel: [ *{html,txt}              ]
    versions                      = ch_versions                      // channel: [ versions.yml             ]
}
