//
// Run VEP to annotate VCF files
//

include { BCFTOOLS_SPLIT } from '../../../modules/nf-core/bcftools/split/main'
include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'
include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    ch_vcf                      // channel: [ val(meta), path(vcf), [path(custom_file1), path(custom_file2)... (optionnal)]]
    ch_fasta                    // channel: [ val(meta2), path(fasta) ] (optional)
    val_genome                  //   value: genome to use
    val_species                 //   value: species to use
    val_cache_version           //   value: cache version to use
    ch_cache                    // channel: [ val(meta3), path(cache) ] (optional)
    ch_extra_files              // channel: [ path(file1), path(file2)... ] (optional)

    main:
    ch_versions = Channel.empty()

    ENSEMBLVEP_VEP(
        ch_vcf,
        val_genome,
        val_species,
        val_cache_version,
        ch_cache,
        ch_fasta,
        ch_extra_files
    )

    TABIX_TABIX(ENSEMBLVEP_VEP.out.vcf)

    ENSEMBLVEP_VEP(
        ch_split_vcf,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        fasta,
        vep_extra_files
    )
    ch_split_annotated_vcf = ENSEMBLVEP_VEP.out.vcf
        .map{meta, vcf ->
            return [
                groupKey(
                    meta - meta.subMap('chunks'),
                    meta.chunks
                ),
                vcf
            ]
        }
        .groupTuple(by: [0])
        .map{meta, vcf ->
            return [meta, vcf.flatten(), []}

    BCFTOOLS_MERGE(ch_split_annotated_vcf)

    TABIX_TABIX(BCFTOOLS_MERGE.out.merged_variants)
    ch_vcf_tbi = BCFTOOLS_MERGE.out.merged_variants
        .join(
            TABIX_TABIX.out.tbi,
            failOnDuplicate: true,
            failOnMismatch: true
        )

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BCFTOOLS_SPLIT.out.versions)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), path(vcf), path(tbi) ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), path(json) ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), path(tab) ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ path(html) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
