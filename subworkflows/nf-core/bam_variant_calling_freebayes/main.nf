include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { FREEBAYES      } from '../../../modules/nf-core/freebayes/main'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_SORT  } from '../../../modules/nf-core/bcftools/sort/main'

workflow BAM_VARIANT_CALLING_FREEBAYES {

    take:
    ch_input        // channel: [ val(meta), bam/cram, bai/crai, bam/cram, bai/crai, bed ]
    genome          // channel: /path/to/fasta
    samples         // channel: /path/to/samples
    populations     // channel: /path/to/populations
    cnv             // channel: /path/to/cnv

    main:

    ch_versions = Channel.empty()

    // Generate index for genome fasta
    ch_genome = Channel.fromPath ( genome ).map { file -> [ [ id: file.baseName ], file ] }
    SAMTOOLS_FAIDX ( ch_genome )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    // Variant calling
    FREEBAYES ( ch_input, genome, SAMTOOLS_FAIDX.out.fai.map { meta, file -> file }, samples, populations, cnv )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())

    // Sort VCF files
    BCFTOOLS_SORT ( FREEBAYES.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    // Index VCF files
    BCFTOOLS_INDEX ( BCFTOOLS_SORT.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), vcf ]
    csi      = BCFTOOLS_INDEX.out.csi          // channel: [ val(meta), csi ]
    tbi      = BCFTOOLS_INDEX.out.tbi          // channel: [ val(meta), tbi ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

