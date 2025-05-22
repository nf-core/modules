include { FREEBAYES      } from '../../../modules/nf-core/freebayes'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_SORT  } from '../../../modules/nf-core/bcftools/sort'

workflow BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS {
    take:
    ch_input       // channel: [mandatory] [ val(meta), path(input1), path(index1), path(input2), path(index2), path(bed) ]
    ch_fasta_fai   // channel: [mandatory] [ val(meta2), path(fasta), path(fai) ]
    ch_samples     // channel: [optional]  [ val(meta3), path(samples) ]
    ch_populations // channel: [optional]  [ val(meta4), path(populations) ]
    ch_cnv         // channel: [optional]  [ val(meta5), path(cnv) ]

    main:

    versions = Channel.empty()

    // Variant calling
    FREEBAYES(
        ch_input,
        ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        ch_fasta_fai.map { meta, _fasta, fai -> [meta, fai] },
        ch_samples,
        ch_populations,
        ch_cnv,
    )

    // Sort VCF files
    BCFTOOLS_SORT(FREEBAYES.out.vcf)

    // Index VCF files
    BCFTOOLS_INDEX(BCFTOOLS_SORT.out.vcf)

    versions = versions.mix(BCFTOOLS_INDEX.out.versions)
    versions = versions.mix(BCFTOOLS_SORT.out.versions)
    versions = versions.mix(FREEBAYES.out.versions)

    emit:
    csi      = BCFTOOLS_INDEX.out.csi // channel: [ val(meta), path(csi) ]
    tbi      = BCFTOOLS_INDEX.out.tbi // channel: [ val(meta), path(tbi) ]
    vcf      = BCFTOOLS_SORT.out.vcf  // channel: [ val(meta), path(vcf) ]
    versions                          // channel: [ path(versions.yml) ]
}
