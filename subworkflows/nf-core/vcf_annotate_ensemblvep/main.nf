//
// Run VEP to annotate VCF files
//

include { BCFTOOLS_SPLIT } from '../../../modules/nf-core/bcftools/split/main'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
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

    BCFTOOLS_SPLIT(ch_vcf)
    ch_split_vcf = BCFTOOLS_SPLIT.out.split_vcf
        .map{ meta, vcf -> [ meta + [ size:vcf instanceof List ? vcf.size() : 1 ], vcf instanceof List ? vcf : [ vcf ] ] }.transpose()

    ENSEMBLVEP_VEP(
        ch_split_vcf,
        val_vep_genome,
        val_vep_species,
        val_vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files
    )
    ch_split_annotated_vcf = ENSEMBLVEP_VEP.out.vcf
        .map{meta, vcf ->
            return [groupKey(meta - meta.subMap('chunks'), meta.chunks), vcf]
        }
        .groupTuple(by: [0])
        .branch{ meta, vcf ->
            single: vcf !instanceof List || vcf instanceof List && vcf.size() == 1
                return [meta, [vcf].flatten()]
            multi:  vcf instanceof List && vcf.size() > 1
                return [meta, vcf.flatten(), []]
        }

    BCFTOOLS_CONCAT(ch_split_annotated_vcf.multi)

    ch_merged_split_vcf = ch_split_annotated_vcf.single
        .mix(BCFTOOLS_CONCAT.out.vcf)

    TABIX_TABIX(ch_merged_split_vcf)

    ch_vcf_tbi = ch_merged_split_vcf
        .join(
            TABIX_TABIX.out.tbi,
            failOnDuplicate: true,
            failOnMismatch: true
        )

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BCFTOOLS_SPLIT.out.versions)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), json ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), tab ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ *.html ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
