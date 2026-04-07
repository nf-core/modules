process SCVITOOLS_SOLO {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.3-' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8764e4208e9639a54d636fc65c839c55dedbfd68def57baea90d1d2007d6a7f/data':
        'community.wave.seqera.io/library/scvi-tools:1.3.3--df115aabdccb7d6b' }"

    input:
    tuple val(meta), path(h5ad)
    val(batch_key)
    val(max_epochs)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.pkl") , emit: predictions
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'solo.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=./tmp

    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scvi: \$(python3 -c 'import scvi; print(scvi.__version__)')
    END_VERSIONS
    """
}
