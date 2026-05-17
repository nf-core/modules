process VIRUSRECOM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/virusrecom:1.4.0--pyhdfd78af_0':
        'quay.io/biocontainers/virusrecom:1.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(alignment), path(mapping)

    output:
    tuple val(meta), path("${prefix}/"), emit: results
    tuple val("${task.process}"), val('virusrecom'), eval("virusrecom --version 2>&1 | head -1"), topic: versions, emit: versions_virusrecom

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    virusrecom \\
        -a ${alignment} \\
        -map ${mapping} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/stub_output.txt
    """
}
