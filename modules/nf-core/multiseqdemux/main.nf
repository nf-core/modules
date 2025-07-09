process MULTISEQDEMUX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/r-seurat_r-seuratobject:4c5a804804327d29'
        : 'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827'}"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
    tuple val(meta), path("*_params_multiseqdemux.csv"), emit: params
    tuple val(meta), path("*_res_multiseqdemux.csv")   , emit: results
    tuple val(meta), path("*_multiseqdemux.rds")       , emit: rds
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('MultiSeqDemux.R')

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_params_multiseqdemux.csv
    touch ${prefix}_res_multiseqdemux.csv
    touch ${prefix}_multiseqdemux.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
    END_VERSIONS
    """
}
