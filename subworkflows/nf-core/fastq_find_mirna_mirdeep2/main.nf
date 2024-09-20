include { BOWTIE_BUILD                               } from '../../../modules/nf-core/bowtie/build/main'
include { SEQKIT_FQ2FA                               } from '../../../modules/nf-core/seqkit/fq2fa/main'
include { SEQKIT_REPLACE                             } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_GENOME    } from '../../../modules/nf-core/seqkit/replace/main'
include { MIRDEEP2_MAPPER                            } from '../../../modules/nf-core/mirdeep2/mapper/main'
include { MIRDEEP2_MIRDEEP2                          } from '../../../modules/nf-core/mirdeep2/mirdeep2/main'

workflow FASTQ_FIND_MIRNA_MIRDEEP2 {

    take:
    ch_reads                     // channel: [ val(meta),  fastq  ]
    ch_genome_fasta              // channel: [ val(meta),  genome_fasta  ]
    ch_mirna_mature_hairpin      // channel: [ val(meta),  mature_mirna, hairpin_mirna ]

    main:

    ch_versions = Channel.empty()

    SEQKIT_FQ2FA ( ch_reads )
    ch_versions = ch_versions.mix(SEQKIT_FQ2FA.out.versions)

    SEQKIT_REPLACE ( SEQKIT_FQ2FA.out.fasta )
    ch_versions = ch_versions.mix(SEQKIT_REPLACE.out.versions)

    SEQKIT_REPLACE_GENOME ( ch_genome_fasta )
    ch_versions = ch_versions.mix(SEQKIT_REPLACE_GENOME.out.versions)

    BOWTIE_BUILD ( SEQKIT_REPLACE_GENOME.out.fastx )
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)

    MIRDEEP2_MAPPER ( SEQKIT_REPLACE.out.fastx, BOWTIE_BUILD.out.index )
    ch_versions = ch_versions.mix(MIRDEEP2_MAPPER.out.versions)

    MIRDEEP2_MIRDEEP2 ( MIRDEEP2_MAPPER.out.mirdeep2_inputs, SEQKIT_REPLACE_GENOME.out.fastx, ch_mirna_mature_hairpin )
    ch_versions = ch_versions.mix(MIRDEEP2_MIRDEEP2.out.versions)

    emit:
    outputs        = MIRDEEP2_MIRDEEP2.out.outputs         // channel: [ val(meta), [ bed, csv, html ] ]
    versions       = ch_versions                           // channel: [ versions.yml ]
}