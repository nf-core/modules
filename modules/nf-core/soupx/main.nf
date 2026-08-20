process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7934a95be2b2e704032cadba8f19685f32d1cebc8febad71e20b43c4f896a7f/data':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-base_r-leidenbase_pruned:4e3c0f41d63a217a' }"

    input:
    tuple val(meta), path(h5ad_filtered), path(h5ad_raw)
    val(npcs)
    val(cluster_algorithm)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions_soupx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_soupx"
    if ("${h5ad_filtered}" == "${prefix}.h5ad" || "${h5ad_raw}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template 'soupx.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_soupx"
    if ("${h5ad_filtered}" == "${prefix}.h5ad" || "${h5ad_raw}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        soupx: \$(Rscript -e "suppressMessages(library(SoupX)); cat(as.character(packageVersion('SoupX')))")
        anndatar: \$(Rscript -e "suppressMessages(library(anndataR)); cat(as.character(packageVersion('anndataR')))")
        seurat: \$(Rscript -e "suppressMessages(library(Seurat)); cat(as.character(packageVersion('Seurat')))")
        leidenbase: \$(Rscript -e "suppressMessages(library(leidenbase)); cat(as.character(packageVersion('leidenbase')))")
    END_VERSIONS
    """
}
