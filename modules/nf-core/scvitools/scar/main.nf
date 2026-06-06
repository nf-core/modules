process SCVITOOLS_SCAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'ghcr.io/scverse/scvi-tools:py3.13-cu12-1.4.3-' :
        workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/df/dfb4b54fd5cb5c5624d947914f5ab8ac86eeaa5f762ebc52c034fbd36cf30250/data':
        'community.wave.seqera.io/library/scvi-tools:1.4.3--cce8c95b58ececa6' }"

    input:
    tuple val(meta), path(filtered), path(unfiltered)
    val(input_layer)
    val(output_layer)
    val(max_epochs)
    val(n_batch)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

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
