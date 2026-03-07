include { SEQFU_STATS as SEQFU_STATS_AFTER  } from '../../../modules/nf-core/seqfu/stats'
include { SEQFU_STATS as SEQFU_STATS_BEFORE } from '../../../modules/nf-core/seqfu/stats'
include { SEQKIT_REPLACE                    } from '../../../modules/nf-core/seqkit/replace'
include { SEQKIT_RMDUP                      } from '../../../modules/nf-core/seqkit/rmdup'
include { SEQKIT_SEQ                        } from '../../../modules/nf-core/seqkit/seq'

workflow FAA_SEQFU_SEQKIT {
    take:
    ch_fasta // tuple val(meta), path(fasta)
    skip_preprocessing // boolean

    main:
    ch_multiqc_files = channel.empty()
    ch_versions = channel.empty()

    SEQFU_STATS_BEFORE(ch_fasta)
    ch_multiqc_files = ch_multiqc_files.mix(SEQFU_STATS_BEFORE.out.multiqc)

    if (!skip_preprocessing) {
        SEQKIT_SEQ(ch_fasta)
        ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

        SEQKIT_REPLACE(SEQKIT_SEQ.out.fastx)

        SEQKIT_RMDUP(SEQKIT_REPLACE.out.fastx)
        ch_fasta = SEQKIT_RMDUP.out.fastx
        ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions)

        SEQFU_STATS_AFTER(SEQKIT_RMDUP.out.fastx)
        ch_multiqc_files = ch_multiqc_files.mix(SEQFU_STATS_AFTER.out.multiqc)
    }

    emit:
    fasta         = ch_fasta
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
