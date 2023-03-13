//
// Run VEP to annotate VCF files
//

include { BCFTOOLS_SPLIT } from '../../../modules/nf-core/bcftools/split/main'
include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'
include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    vcf               // channel: [ val(meta), vcf, tbi ]
    fasta             //   value: fasta to use (optionnal)
    vep_genome        //   value: genome to use
    vep_species       //   value: species to use
    vep_cache_version //   value: cache version to use
    vep_cache         //    path: /path/to/vep/cache (optionnal)
    vep_extra_files   // channel: [ file1, file2...] (optionnal)

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_SPLIT(vcf)
    ch_split_vcf = BCFTOOLS_SPLIT.out.vcf
        .map{meta, vcf ->
            meta.chunks = vcf instanceof List ? vcf.size() : 1
            return [meta, vcf]
        }
        .transpose(by: [0])

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
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), json ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), tab ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ *.html ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
