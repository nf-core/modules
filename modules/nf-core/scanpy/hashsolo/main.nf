process SCANPY_HASHSOLO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/168ecbbe27ccef766741ccbf937b0d2675be2e19b0565035e0719f1e9ea5ee95/data':
        'community.wave.seqera.io/library/python_pyyaml_scanpy:b5509a698e9aae25' }"

    input:
    tuple val(meta), path(input_h5ad), val(cell_hashing_columns)
    val(priors)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('hashsolo.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    # Prevent failures during scanpy import
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -f 2 -d " ")
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
