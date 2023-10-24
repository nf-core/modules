include { FREEBAYES      } from '../../../modules/nf-core/freebayes/main'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_SORT  } from '../../../modules/nf-core/bcftools/sort/main'

workflow BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS {

    take:
    ch_input        // channel: [mandatory] [ val(meta), path(input1), path(index1), path(input2), path(index2), path(bed) ]
    ch_fasta_fai    // channel: [mandatory] [ val(meta2), path(fasta), path(fai) ]
    ch_samples      // channel: [optional]  [ path(samples) ]
    ch_populations  // channel: [optional]  [ path(populations ]
    ch_cnv          // channel: [optional]  [ path(cnv) ]

    main:

    ch_versions = Channel.empty()

    // Variant calling
    FREEBAYES ( ch_input, ch_fasta_fai.map{ meta, fasta, fai -> fasta }, ch_fasta_fai.map{ meta, fasta, fai -> fai }, ch_samples, ch_populations, ch_cnv )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())

    // Sort VCF files
    BCFTOOLS_SORT ( FREEBAYES.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    // Index VCF files
    BCFTOOLS_INDEX ( BCFTOOLS_SORT.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), path(vcf) ]
    csi      = BCFTOOLS_INDEX.out.csi          // channel: [ val(meta), path(csi) ]
    tbi      = BCFTOOLS_INDEX.out.tbi          // channel: [ val(meta), path(tbi) ]

    versions = ch_versions                     // channel: [ path(versions.yml) ]
}

