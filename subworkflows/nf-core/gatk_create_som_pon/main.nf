//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

params.mutect2_options      = [:]
params.gendbimport_options  = [:]
params.createsompon_options = [:]

include { GATK4_MUTECT2                     } from '../../../modules/gatk/mutect2/main'                     addParams( options: params.mutect2_options )
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/gatk/genomicsdbimport/main'            addParams( options: params.gendbimport_options )
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/gatk/createsomaticpanelofnormals/main' addParams( options: params.createsompon_options )

workflow GATK_CREATE_SOM_PON {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai/csi] ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_STATS ( ch_bam_bai, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
