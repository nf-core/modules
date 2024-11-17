process DOUBLETDETECTION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b7a5183db20eee4e6eed6b2ec6db7202023b484debd48abcd273cb4658cdc80c/data' :
        'community.wave.seqera.io/library/anndata_louvain_pip_doubletdetection:cbe92394c10372fa' }"

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
