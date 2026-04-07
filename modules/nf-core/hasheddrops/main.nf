process HASHEDDROPS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/431b27926fac88d6334ee3e8f63479f69a1a69340b305a05b70bc84083d301aa/data':
        'community.wave.seqera.io/library/bioconductor-dropletutils_r-seurat:e1dff3a0fb7c5920' }"

    input:
    tuple val(meta), path(hto_matrix), val(runEmptyDrops), path(rna_matrix)

    output:
    tuple val(meta), path("*_emptyDrops.png")         , emit: empty_drops_plot
    tuple val(meta), path("*_emptyDrops.csv")         , emit: empty_drops_csv
    tuple val(meta), path("*_emptyDrops.rds")         , emit: empty_drops_rds
    tuple val(meta), path("*_results_hasheddrops.csv"), emit: results
    tuple val(meta), path("*_id_to_hash.csv")         , emit: id_to_hash
    tuple val(meta), path("*_hasheddrops.rds")        , emit: rds
    tuple val(meta), path("*_plot_hasheddrops.png")   , emit: plot
    tuple val(meta), path("*_params_hasheddrops.csv") , emit: params
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('HashedDrops.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_emptyDrops.png
    touch ${prefix}_emptyDrops.csv
    touch ${prefix}_emptyDrops.rds
    touch ${prefix}_results_hasheddrops.csv
    touch ${prefix}_id_to_hash.csv
    touch ${prefix}_hasheddrops.rds
    touch ${prefix}_plot_hasheddrops.png
    touch ${prefix}_params_hasheddrops.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
        dropletutils: \$(Rscript -e "library(DropletUtils); cat(as.character(packageVersion('DropletUtils')))")
    END_VERSIONS
    """
}
