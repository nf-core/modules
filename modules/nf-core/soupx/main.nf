process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-leidenbase_r-seurat_r-soupx:b67b9fe9c42610c9':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-leidenbase_r-seurat_r-soupx:9582ccf4a756d29e' }"

    input:
    tuple val(meta), path(h5ad_filtered), path(h5ad_raw)

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
    def args_map = [:]
    (task.ext.args ?: '').findAll(/--(\S+)\s+(\S+)/) { full, k, v -> args_map[k] = v }
    npcs              = (args_map.npcs ?: 50) as Integer
    cluster_algorithm = [ louvain: 1, louvain_multilevel: 2, slm: 3, leiden: 4 ][args_map.cluster_algorithm ?: 'leiden'] ?: error("Unknown cluster_algorithm '${args_map.cluster_algorithm}'. Must be one of: louvain, louvain_multilevel, slm, leiden.")
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
