process VIRUSRECOM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/virusrecom:1.4.0--pyhdfd78af_0':
        'quay.io/biocontainers/virusrecom:1.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(alignment), path(mapping), path(iwic)
    val(query)

    output:
    tuple val(meta), path("${meta.id}/"), emit: results
    tuple val("${task.process}"), val('virusrecom'), eval("virusrecom -h 2>&1 | sed -n 's/.*Version: \\([^ ]*\\).*/\\1/p'"), topic: versions, emit: versions_virusrecom

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def use_iwic     = iwic ? true : false
    def input_cmd    = use_iwic ? "-iwic ${iwic}" : "-a ${alignment}"
    def map_cmd      = use_iwic ? "" : "-map ${mapping}"
    """
    virusrecom \\
        ${input_cmd} \\
        ${map_cmd} \\
        -q ${query} \\
        -o ${meta.id} \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/stub_output.txt
    """
}
