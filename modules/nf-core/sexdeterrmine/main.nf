process SEXDETERRMINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sexdeterrmine:1.1.2--hdfd78af_1'
        : 'quay.io/biocontainers/sexdeterrmine:1.1.2--hdfd78af_1'}"

    input:
    tuple val(meta), path(depth)
    path sample_list_file

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('sexdeterrmine'), eval("sexdeterrmine --version 2>&1"), topic: versions, emit: versions_sexdeterrmine

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sexdeterrmine"
    def sample_list = sample_list_file ? '-f ${sample_list_file}' : ''
    if ("${depth}" == "${prefix}.tsv") {
        error("Input depth and output TSV names are the same, set prefix in module configuration to disambiguate!")
    }

    """
    sexdeterrmine \\
        -I ${depth} \\
        ${sample_list} \\
        ${args} \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_sexdeterrmine"
    if ("${depth}" == "${prefix}.tsv") {
        error("Input depth and output TSV names are the same, set prefix in module configuration to disambiguate!")
    }

    """
    touch ${prefix}.tsv
    touch ${prefix}.json
    """
}
