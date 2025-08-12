process SCANPY_HASHSOLO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a06fe0be02e1a82be29f415a554cfd3ad9d921cd63e71ebce3b141a723564426/data':
        'community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scanpy:bf7b7ef27120d15c' }"

    input:
    tuple val(meta), path(data), val(cell_hashing_columns)

    output:
    tuple val(meta), path("*_assignment_hashsolo.csv"), emit: assignment
    tuple val(meta), path("*_hashsolo.h5ad")          , emit: h5ad
    tuple val(meta), path("*_params_hashsolo.csv")    , emit: params
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    template('hashsolo.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_assignment_hashsolo.csv
    touch ${prefix}_hashsolo.h5ad
    touch ${prefix}_params_hashsolo.csv

    # Prevent failures during scanpy import
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -f 2 -d " ")
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
