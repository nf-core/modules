process SYLPHTAX_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph-tax:1.9.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/sylph-tax:1.9.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(sylphtax_reports)
    val data_type

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('sylph-tax'), eval("sylph-tax --version 2>&1 | tail -1"), emit: versions_sylphtax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    sylph-tax \\
        merge \\
        ${sylphtax_reports} \\
        --column ${data_type} \\
        --output ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    touch ${prefix}.tsv
    """
}
