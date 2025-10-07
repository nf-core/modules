include { SEQFU_STATS as SEQFU_STATS_BEFORE } from '../../../modules/nf-core/seqfu/stats/main'
include { SEQKIT_SEQ                        } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_RMDUP                      } from '../../../modules/nf-core/seqkit/rmdup/main'
include { SEQKIT_REPLACE                    } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQFU_STATS as SEQFU_STATS_AFTER  } from '../../../modules/nf-core/seqfu/stats/main'

workflow FAA_SEQFU_SEQKIT {

    take:
    fasta              // tuple val(meta), path(fasta)
    skip_preprocessing // boolean

    main:
    ch_fasta         = fasta
    ch_multiqc_files = Channel.empty()
    ch_versions      = Channel.empty()

    SEQFU_STATS_BEFORE( ch_fasta )
    ch_multiqc_files = ch_multiqc_files.mix( SEQFU_STATS_BEFORE.out.multiqc )
    ch_versions      = ch_versions.mix( SEQFU_STATS_BEFORE.out.versions )

    if (!skip_preprocessing) {
        SEQKIT_SEQ( ch_fasta )
        ch_versions = ch_versions.mix( SEQKIT_SEQ.out.versions )

        SEQKIT_REPLACE( SEQKIT_SEQ.out.fastx )
        ch_versions = ch_versions.mix( SEQKIT_REPLACE.out.versions )

        SEQKIT_RMDUP( SEQKIT_REPLACE.out.fastx )
        ch_fasta    = SEQKIT_RMDUP.out.fastx
        ch_versions = ch_versions.mix( SEQKIT_RMDUP.out.versions )

        SEQFU_STATS_AFTER( SEQKIT_RMDUP.out.fastx )
        ch_multiqc_files = ch_multiqc_files.mix( SEQFU_STATS_AFTER.out.multiqc )
        ch_versions      = ch_versions.mix( SEQFU_STATS_AFTER.out.versions )
    }

    emit:
    fasta         = ch_fasta
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
