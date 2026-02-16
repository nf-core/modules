//
// Run VEP and/or SNPEFF to annotate VCF files
//

include { BCFTOOLS_CONCAT        } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_PLUGINSCATTER } from '../../../modules/nf-core/bcftools/pluginscatter'
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort'
include { ENSEMBLVEP_VEP         } from '../../../modules/nf-core/ensemblvep/vep'
include { SNPEFF_SNPEFF          } from '../../../modules/nf-core/snpeff/snpeff'
include { TABIX_BGZIP            } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_TABIX            } from '../../../modules/nf-core/tabix/tabix'

workflow VCF_ANNOTATE_ENSEMBLVEP_SNPEFF {
    take:
    ch_vcf // channel: [ val(meta), path(vcf), path(tbi), [path(file1), path(file2)...] ]
    ch_fasta // channel: [ val(meta2), path(fasta) ] (optional)
    val_vep_genome //   value: genome to use
    val_vep_species //   value: species to use
    val_vep_cache_version //   value: cache version to use
    ch_vep_cache // channel: [ path(cache) ] (optional)
    ch_vep_extra_files // channel: [ path(file1), path(file2)... ] (optional)
    val_snpeff_db //   value: the db version to use for snpEff
    ch_snpeff_cache // channel: [ path(cache) ] (optional)
    val_tools_to_use //   value: a list of tools to use options are: ["ensemblvep", "snpeff"]
    val_sites_per_chunk //   value: the amount of variants per scattered VCF

    main:
    def ch_versions = channel.empty()
    def ch_vep_input = channel.empty()
    def ch_scatter = channel.empty()

    // Scatter if val_sites_per_chunk is set up
    if (val_sites_per_chunk) {
        // Prepare the input VCF channel for scattering (split VCFs from custom files)
        def ch_input = ch_vcf.multiMap { meta, vcf, tbi, custom_files ->
            vcf: [meta, vcf, tbi]
            custom: [meta, custom_files]
        }

        // Scatter the input VCFs into multiple VCFs
        // These VCFs contain the amount of variants specified by `val_sites_per_chunk`
        // The lower this value is, the more files will be created

        BCFTOOLS_PLUGINSCATTER(
            ch_input.vcf,
            val_sites_per_chunk,
            [],
            [],
            [],
            [],
        )

        // If BCFTOOLS_PLUGINSCATTER created multiple files we return a list of vcfs and the size of that list
        // Otherwise, a single vcf and the value 1
        // We then use transpose and combine to
        // - Transpose on the VCFs => Creates an entry for each VCF in the list
        // - Re-add the sample specific custom files
        // multiMap is used to
        // - Define the new ID. The `_annotated` suffix is used to disambiguate the VEP output with its input
        // - Create 2 channels: one with the VEP input and one with the original ID and count of scattered VCFs
        ch_scatter = BCFTOOLS_PLUGINSCATTER.out.scatter
            .map { meta, vcfs ->
                def is_list = vcfs instanceof ArrayList
                count = is_list ? vcfs.size() : 1
                [meta, is_list ? vcfs : [vcfs], count]
            }
            .transpose(by: 1)
            .combine(ch_input.custom, by: 0)
            .multiMap { meta, vcf, count, custom_files ->
                def new_id = "${meta.id}${vcf.name.replace(meta.id, "").tokenize(".")[0]}_annotated" as String
                input: [meta + [id: new_id], vcf, custom_files]
                count: [meta + [id: new_id], meta.id, count]
            }

        ch_vep_input = ch_scatter.input
    }
    else {
        // Use the normal input when no scattering has to be performed
        ch_vep_input = ch_vcf.map { meta, vcf, _tbi, files -> [meta, vcf, files] }
    }

    // Annotate with ensemblvep if it's part of the requested tools
    def ch_vep_output = channel.empty()
    def ch_vep_reports = channel.empty()

    if ("ensemblvep" in val_tools_to_use) {
        ENSEMBLVEP_VEP(
            ch_vep_input,
            val_vep_genome,
            val_vep_species,
            val_vep_cache_version,
            ch_vep_cache,
            ch_fasta,
            ch_vep_extra_files,
        )
        ch_vep_output = ENSEMBLVEP_VEP.out.vcf
        ch_vep_reports = ENSEMBLVEP_VEP.out.report
    }
    else {
        ch_vep_output = ch_vep_input.map { meta, vcf, _files -> [meta, vcf] }
    }

    // Annotate with snpeff if it's part of the requested tools
    def ch_snpeff_output = channel.empty()
    def ch_snpeff_reports = channel.empty()
    def ch_snpeff_html = channel.empty()
    def ch_snpeff_genes = channel.empty()

    if ("snpeff" in val_tools_to_use) {
        SNPEFF_SNPEFF(
            ch_vep_output,
            val_snpeff_db,
            ch_snpeff_cache,
        )
        ch_snpeff_reports = SNPEFF_SNPEFF.out.report
        ch_snpeff_html = SNPEFF_SNPEFF.out.summary_html
        ch_snpeff_genes = SNPEFF_SNPEFF.out.genes_txt

        TABIX_BGZIP(
            SNPEFF_SNPEFF.out.vcf
        )

        ch_snpeff_output = TABIX_BGZIP.out.output
    }
    else {
        ch_snpeff_output = ch_vep_output
    }

    // Gather the files back together if they were scattered
    def ch_ready_vcfs = channel.empty()

    if (val_sites_per_chunk) {
        // Concatenate the VCFs back together with bcftools concat

        def ch_concat_input = ch_snpeff_output
            .join(ch_scatter.count, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, id, count ->
                def new_meta = meta + [id: id]
                [groupKey(new_meta, count), vcf]
            }
            .groupTuple()
            .map { meta, vcf ->
                [meta, vcf, []]
            }

        BCFTOOLS_CONCAT(ch_concat_input)

        // Sort the concatenate output (bcftools concat is unable to do this on its own)

        BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)

        ch_ready_vcfs = BCFTOOLS_SORT.out.vcf
    }
    else {
        ch_ready_vcfs = ch_snpeff_output
    }

    //
    // Index the resulting bgzipped VCFs
    //

    // Split the bgzipped VCFs from the unzipped VCFs (only bgzipped VCFs should be indexed)
    def ch_tabix_input = ch_ready_vcfs.branch { meta, vcf ->
        bgzip: vcf.extension == "gz"
        unzip: true
        return [meta, vcf, []]
    }

    TABIX_TABIX(ch_tabix_input.bgzip)

    def ch_vcf_tbi = ch_tabix_input.bgzip
        .join(TABIX_TABIX.out.index, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_tabix_input.unzip)

    emit:
    snpeff_genes   = ch_snpeff_genes // channel: [ val(meta), path(genes) ]
    snpeff_html    = ch_snpeff_html // channel: [ val(meta), path(html) ]
    snpeff_reports = ch_snpeff_reports // channel: [ val(meta), path(csv) ]
    vcf_tbi        = ch_vcf_tbi // channel: [ val(meta), path(vcf), path(tbi) ]
    vep_reports    = ch_vep_reports // channel: [ path(html) ]
    versions       = ch_versions // channel: [ versions.yml ]
}
