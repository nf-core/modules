process SCANPY_SCRUBLET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/876de316db5929c338bf04e65a5b993074acacab9dfe23274fdaa9dafa580081/data'
        : 'community.wave.seqera.io/library/scanpy_scikit-image:7d3d2c63b2063667'}"

    input:
    tuple val(meta), path(h5ad)
    val batch_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.pkl"), emit: predictions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('scrublet.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
