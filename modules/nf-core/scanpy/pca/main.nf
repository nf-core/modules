process SCANPY_PCA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bd4462356828eb975a9ae0a3073ec720d19f4918106b696e3c298a556f288480/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ffeb5af98a9a5bde'}"

    input:
    tuple val(meta), path(h5ad)
    val key_added

    output:
    tuple val(meta), path("*.h5ad") , emit: h5ad
    tuple val(meta), path("X_*.pkl"), emit: obsm
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_pca"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('pca.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_pca"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    # These are needed to prevent errors during import of scanpy
    # when using singularity/apptainer
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    touch ${prefix}.h5ad
    touch X_${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
