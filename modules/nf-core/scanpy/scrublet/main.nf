process SCANPY_SCRUBLET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a37ea9b603007f0751dff6e5d74271c6c268fc2a1bab35bb11af06c2a24d851/data'
        : 'community.wave.seqera.io/library/python_scanpy_scikit-image:9f0ffb2755afe2a1'}"

    input:
    tuple val(meta), path(h5ad)
    val batch_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.pkl") , emit: predictions
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('scrublet.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
