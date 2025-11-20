process SCVITOOLS_SCAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8764e4208e9639a54d636fc65c839c55dedbfd68def57baea90d1d2007d6a7f/data':
        'community.wave.seqera.io/library/scvi-tools:1.3.3--df115aabdccb7d6b' }"

    input:
    tuple val(meta), path(filtered), path(unfiltered)
    val(input_layer)
    val(output_layer)
    val(max_epochs)
    val(n_batch)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'scar.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version)
        scvi: \$(python -c "import scvi; print(scvi.__version__)")
    END_VERSIONS
    """
}
