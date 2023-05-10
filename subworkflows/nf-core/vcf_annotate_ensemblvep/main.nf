//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP         } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_PLUGINSCATTER } from '../../../modules/nf-core/bcftools/pluginscatter/main'
include { BCFTOOLS_CONCAT        } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort/main'

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
    // Prepare the input VCF channel for scattering (split VCFs from custom files)
    //

    ch_input = ch_vcf
        .multiMap { meta, vcf, tbi, custom_files ->
            vcf:    [ meta, vcf, tbi ]
            custom: [ meta, custom_files ]
        }

    //
    // Scatter the input VCFs into multiple VCFs. These VCFs contain the amount of variants
    // specified by `val_sites_per_chunk`. The lower this value is, the more files will be created
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
            is_list = vcfs instanceof ArrayList
            count = is_list ? vcfs.size() : 1
            [ meta, is_list ? vcfs : [vcfs], count ]
            // Channel containing the list of VCFs and the size of this list
        }
        .transpose(by:1) // Transpose on the VCFs => Creates an entry for each VCF in the list
        .combine(ch_input.custom, by: 0) // Re-add the sample specific custom files
        .multiMap { meta, vcf, count, custom_files ->
            // Define the new ID. The `_annotated` is to disambiguate the VEP output with its input
            new_id = "${meta.id}${vcf.name.replace(meta.id,"").tokenize(".")[0]}_annotated" as String
            new_meta = meta + [id:new_id]

            // Create channels: one with the VEP input and one with the original ID and count of scattered VCFs
            vep:    [ new_meta, vcf, custom_files ]
            count:  [ new_meta, meta.id, count ]
        }

    ENSEMBLVEP_VEP(
        ch_scatter.vep,
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
        .join(ch_scatter.count, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, id, count ->
            new_meta = meta + [id:id]
            [ groupKey(new_meta, count), vcf ]
        }
        .groupTuple() // Group the VCFs which need to be concatenated
        .map { it + [[]] }

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    //
    // Sort the concatenate output (bcftools concat is unable to do this on its own)
    //

    BCFTOOLS_SORT(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    //
    // Index the resulting bgzipped VCFs
    //

    ch_tabix_input = BCFTOOLS_SORT.out.vcf
        .branch { meta, vcf ->
            // Split the bgzipped VCFs from the unzipped VCFs (only bgzipped VCFs should be indexed)
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
