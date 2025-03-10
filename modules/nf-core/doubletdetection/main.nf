process DOUBLETDETECTION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/878a0582c1de0ad7370ad1fbdebd7e786c77d29b064e10a7c09c35a9df3bfb97/data' :
        'community.wave.seqera.io/library/anndata_louvain_numpy_pip_pruned:9ff7bfd3c5201947' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.pkl") , emit: predictions
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'doubletdetection.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=./tmp
    export NUMBA_CACHE_DIR=./tmp

    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 -c 'import platform as pf; print(pf.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        doubletdetection: \$(python3 -c 'import doubletdetection as dt; print(dt.__version__)')
    END_VERSIONS
    """
}
