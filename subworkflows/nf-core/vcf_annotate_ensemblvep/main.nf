//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_QUERY            } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_MERGE            } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_PLUGINSCATTER    } from '../../../modules/nf-core/bcftools/pluginscatter/main'
include { BCFTOOLS_CONCAT           } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_PLUGINSPLIT      } from '../../../modules/nf-core/bcftools/pluginsplit/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    ch_vcf                      // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_fasta                    // channel: [ val(meta2), path(fasta) ] (optional)
    ch_fasta_fai                // channel: [ val(meta2), path(fasta_fai) ] (optional)
    val_genome                  //   value: genome to use
    val_species                 //   value: species to use
    val_cache_version           //   value: cache version to use
    ch_cache                    // channel: [ path(cache) ] (optional)
    ch_extra_files              // channel: [ path(file1), path(file2)... ] (optional)
    val_sites_per_chunk         //   value: the amount of variants per scattered VCF

    main:
    ch_versions = Channel.empty()

    //
    // Define the sample names in each VCF using bcftools query --list-samples
    //

    ch_vcf_map = ch_vcf
        .multiMap { meta, vcf, tbi ->
            vcf:   [ meta, vcf, tbi ]
            meta:  [ meta.id, meta ]
        }

    BCFTOOLS_QUERY(
        ch_vcf_map.vcf,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    ch_samples_list = BCFTOOLS_QUERY.out.txt
        .map { meta, txt ->
            return create_samples_file(meta, txt)
        }
        .collectFile(name:"samples.txt", newLine:true)

    //
    // Merge the VCFs into one big VCF using bcftools merge
    //

    ch_merge_input = ch_vcf_map.vcf
        .map { meta, vcf, tbi ->
            [ [id:'raw'], vcf, tbi ]
        }
        .groupTuple() // No size needed here because this will merge all VCFs in the workflow

    BCFTOOLS_MERGE(
        ch_merge_input,
        [],
        ch_fasta instanceof List ? [] : ch_fasta.map { it[1] },
        ch_fasta_fai instanceof List ? [] : ch_fasta_fai.map { it[1] }
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    //
    // Scatter the merged VCF into multiple smaller VCFs
    //

    BCFTOOLS_PLUGINSCATTER(
        BCFTOOLS_MERGE.out.merged_variants.map { it + [[]] },
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

    ch_vep_input = BCFTOOLS_PLUGINSCATTER.out.scatter
        .transpose()
        .map { meta, vcf ->
            [ [id:vcf.simpleName], vcf, [] ]
        }

    ENSEMBLVEP_VEP(
        ch_vep_input,
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
        .map { meta, vcf ->
            [ [id:'annotated'], vcf ]
        }
        .groupTuple() // Again no size needed here because it will group all files in the workflow
        .map { it + [[]] }

    BCFTOOLS_CONCAT(
        ch_concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    //
    // Split the VCFs back into their original forms
    //

    BCFTOOLS_PLUGINSPLIT(
        BCFTOOLS_CONCAT.out.vcf.map { it + [[]] },
        ch_samples_list,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    //
    // Index the resulting gzipped VCFs
    //

    ch_tabix_input = BCFTOOLS_PLUGINSPLIT.out.vcf
        .transpose()
        .map { meta, vcf ->
            name = vcf.name.replace(".bcf", "").replace(".vcf", "").replace(".gz", "")
            [ name, vcf ]
        } // Re-add the original meta to each VCF
        .join(ch_vcf_map.meta, failOnDuplicate:true, failOnMismatch:true)
        .branch { id, vcf, meta ->
            bgzip: vcf.extension == "gz"
                return [ meta, vcf ]
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

def create_samples_file(meta, samples) {
    def file_prefix = meta.id
    def samples_list = samples.readLines().join(",")
    return "${samples_list}\t-\t${file_prefix}"
}
