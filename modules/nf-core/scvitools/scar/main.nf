process SCVITOOLS_SCAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata_scvi-tools:54d2eb2f946e0e43':
        'community.wave.seqera.io/library/anndata_scvi-tools:fa9451a13918eae0' }"

    input:
    tuple val(meta), path(h5ad), path(raw)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'scar.py'

    stub:
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
        scvi-tools: \$(python -c "import scvi; print(scvi.__version__)")
    END_VERSIONS
    """
}
