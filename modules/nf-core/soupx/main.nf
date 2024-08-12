process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_leidenalg_python-igraph_pruned:c2538291aadd50cb':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_leidenalg_python-igraph_pruned:1b1b2ad4205f41be' }"

    input:
    tuple val(meta), path(filtered), path(unfiltered)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'soupx.py'

    stub:
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata2ri: \$(python -c "import anndata2ri; print(anndata2ri.__version__)")
        rpy2: \$(python -c "import rpy2; print(rpy2.__version__)")
        SoupX: \$(Rscript -e "library(soupx); cat(as.character(packageVersion('soupx')))")
    END_VERSIONS
    """
}
