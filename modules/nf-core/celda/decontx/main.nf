process CELDA_DECONTX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-celda_anndata_numpy:2aed5fa978c663d9':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-celda_anndata_numpy:63af229ac9152259' }"

    input:
    tuple val(meta), path(h5ad), path(raw)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    batch_col = task.ext.batch_col ?: "batch"
    template 'decontx.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    batch_col = task.ext.batch_col ?: "batch"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 -c "import platform; print(platform.python_version())")
        anndata: \$(python3 -c "import anndata as ad; print(ad.__version__)")
        anndata2ri: \$(python3 -c "import anndata2ri; print(anndata2ri.__version__)")
        rpy2: \$(python3 -c "import rpy2; print(rpy2.__version__)")
        celda: \$(python3 -c "import anndata2ri; import rpy2; import rpy2.robjects as ro; celda = ro.packages.importr('celda'); print(celda.__version__)")
    END_VERSIONS
    """
}
