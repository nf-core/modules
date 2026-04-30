process BIOMFORMAT_CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biom-format:2.1.15':
        'quay.io/biocontainers/biom-format:2.1.15' }"

    input:
    tuple val(meta), path(biom)

    output:
    tuple val(meta), path("*.biom"), optional: true, emit: biom
    tuple val(meta), path("*.txt") , optional: true, emit: txt
    tuple val("${task.process}"), val('biom-format'), eval("biom --version | sed 's/biom, version //'"), topic: versions, emit: versions_biomformat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = task.ext.args.contains("--to-tsv") ? "${prefix}.txt" : "${prefix}.biom"
    if( "${output}" == "${biom}" ) error "ERROR: Input and output names are the same, set prefix in module configuration"
    """
    biom convert \\
        -i ${biom} \\
        -o ${output} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = task.ext.args.contains("--to-tsv") ? "${prefix}.txt" : "${prefix}.biom"
    if( "${output}" == "${biom}" ) error "ERROR: Input and output names are the same, set prefix in module configuration"
    """
    touch ${output}
    """
}
