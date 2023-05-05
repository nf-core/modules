//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_PLUGINSCATTER    } from '../../../modules/nf-core/bcftools/pluginscatter/main'
include { BCFTOOLS_CONCAT           } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    ch_vcf                      // channel: [ val(meta), path(vcf), path(tbi), [path(file1), path(file2)...] ]
    ch_fasta                    // channel: [ val(meta2), path(fasta) ] (optional)
    val_genome                  //   value: genome to use
    val_species                 //   value: species to use
    val_cache_version           //   value: cache version to use
    ch_cache                    // channel: [ path(cache) ] (optional)
    ch_extra_files              // channel: [ path(file1), path(file2)... ] (optional)
    val_sites_per_chunk         //   value: the amount of variants per scattered VCF

    main:
    ch_versions = Channel.empty()

    //
    // Prepare the input VCF channel for scattering
    //

    ch_input = ch_vcf
        .multiMap { meta, vcf, tbi, custom_files ->
            vcf:    [ meta, vcf, tbi ]
            custom: [ meta, custom_files ]
        }

    //
    // Scatter the input VCFs into multiple VCFs. These VCFs contain the amount of variants
    // specified by `val_sites_per_chunk`.
    //

    BCFTOOLS_PLUGINSCATTER(
        ch_input.vcf,
        val_sites_per_chunk,
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSCATTER.out.versions.first())

    //
    // Run the annotation with EnsemblVEP
    //

    ch_scatter = BCFTOOLS_PLUGINSCATTER.out.scatter
        .map { meta, vcfs ->
            count = vcfs instanceof ArrayList ? vcfs.size() : 1
            [ meta, vcfs instanceof ArrayList ? vcfs : [vcfs], count ] // Channel containing the list of VCFs
        }
        .transpose(by:1) // Transpose on the VCFs
        .combine(ch_input.custom, by: 0)
        .multiMap { meta, vcf, count, custom_files ->
            new_id = "${meta.id}${vcf.name.replace(meta.id,"").tokenize(".")[0]}_annotated" as String
            new_meta = meta + [id:new_id]
            vcf:    [ new_meta, vcf, custom_files ]
            id:     [ new_meta, meta.id ]
            counts: [ new_meta, count ]
        }

    ENSEMBLVEP_VEP(
        ch_scatter.vcf,
        val_genome,
        val_species,
        val_cache_version,
        ch_cache,
        ch_fasta,
        ch_extra_files
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())

    //
    // Concatenate the VCFs back together with bcftools concat
    //

    ch_concat_input = ENSEMBLVEP_VEP.out.vcf
        .join(ch_scatter.id, failOnDuplicate:true, failOnMismatch:true)
        .join(ch_scatter.counts, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, id, count ->
            new_meta = meta + [id:id]
            [ groupKey(new_meta, count), vcf ]
        }
        .groupTuple()
        .map { it + [[]] }

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    BCFTOOLS_SORT(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    //
    // Index the resulting bgzipped VCFs
    //

    ch_tabix_input = BCFTOOLS_SORT.out.vcf
        .branch { meta, vcf ->
            bgzip: vcf.extension == "gz"
            unzip: true
                return [ meta, vcf, [] ]
        }

    TABIX_TABIX(
        ch_tabix_input.bgzip
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    ch_vcf_tbi = ch_tabix_input.bgzip
        .join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_tabix_input.unzip)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), path(vcf), path(tbi) ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ path(html) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
