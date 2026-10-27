include { METHYLDACKEL_EXTRACT } from '../../../modules/nf-core/methyldackel/extract/main'
include { METHYLDACKEL_MBIAS   } from '../../../modules/nf-core/methyldackel/mbias/main'

workflow BAM_METHYLDACKEL {

    take:
    ch_bam   // channel: [ val(meta), bam, bai ]
    ch_fasta // channel: [ val(meta), fasta, fai ]

    main:
    ch_methyldackel_extract_bedgraph  = channel.empty()
    ch_methyldackel_extract_methylkit = channel.empty()
    ch_methyldackel_mbias             = channel.empty()
    ch_multiqc_files                  = channel.empty()

    /*
     * Extract per-base methylation and plot methylation bias
     */

    METHYLDACKEL_EXTRACT (
        ch_bam,
        ch_fasta
    )
    ch_methyldackel_extract_bedgraph  = METHYLDACKEL_EXTRACT.out.bedgraph
    ch_methyldackel_extract_methylkit = METHYLDACKEL_EXTRACT.out.methylkit

    METHYLDACKEL_MBIAS (
        ch_bam,
        ch_fasta
    )
    ch_methyldackel_mbias = METHYLDACKEL_MBIAS.out.txt

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_methyldackel_extract_bedgraph.collect{ _meta, bedgraph -> bedgraph  }
                        .mix(ch_methyldackel_extract_methylkit.collect{ _meta, methylkit -> methylkit })
                        .mix(ch_methyldackel_mbias.collect{ _meta, txt -> txt  })

    emit:
    methyldackel_extract_bedgraph  = ch_methyldackel_extract_bedgraph  // channel: [ val(meta), [ bedgraph ]  ]
    methyldackel_extract_methylkit = ch_methyldackel_extract_methylkit // channel: [ val(meta), [ methylkit ] ]
    methyldackel_mbias             = ch_methyldackel_mbias             // channel: [ val(meta), [ mbias ]     ]
    multiqc                        = ch_multiqc_files                  // channel: [ *{html,txt}              ]
}
