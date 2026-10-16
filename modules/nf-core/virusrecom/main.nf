process VIRUSRECOM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc4e9e01ce104f8f1ba09daf74ec67846a1665a057418d6e06721e9396cf66ea/data'
:         'community.wave.seqera.io/library/virusrecom:1.4.0--b0b2c64800fa6452' }"

    input:
    tuple val(meta), path(alignment), path(mapping), path(iwic)
    val(query)

    output:
    tuple val(meta), path("${prefix}/"), emit: results
    tuple val("${task.process}"), val('virusrecom'), eval("virusrecom -h 2>&1 | sed -n 's/.*Version: \\([^ ]*\\).*/\\1/p'"), topic: versions, emit: versions_virusrecom

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def use_iwic     = iwic ? true : false
    def input_cmd    = use_iwic ? "-iwic ${iwic}" : "-a ${alignment}"
    def map_cmd      = use_iwic ? "" : "-map ${mapping}"
    """
    virusrecom \\
        ${input_cmd} \\
        ${map_cmd} \\
        -q ${query} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/stub_output.txt
    """
}
