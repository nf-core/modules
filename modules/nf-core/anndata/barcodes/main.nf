process ANNDATA_BARCODES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa01776f530c0a0fe2b7d8d41009884040f75d175ff194641f846b215a1da8d6/data':
        'community.wave.seqera.io/library/anndata_pandas_python_pyyaml:3f88bccf1cdade40' }"

    input:
    tuple val(meta), path(h5ad), path(barcodes)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'barcodes.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
    END_VERSIONS
    """
}
