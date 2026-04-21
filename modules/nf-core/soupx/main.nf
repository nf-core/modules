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
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'soupx.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
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
