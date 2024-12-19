process SCVITOOLS_SCAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scvi-tools:1.2.0--adbeb4160de8b08e':
        'community.wave.seqera.io/library/scvi-tools:1.2.0--680d378b86801b8a' }"

    input:
    tuple val(meta), path(filtered), path(unfiltered)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    input_layer = task.ext.input_layer ?: "X"
    output_layer = task.ext.output_layer ?: "scar"
    max_epochs = task.ext.max_epochs ?: ""
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
