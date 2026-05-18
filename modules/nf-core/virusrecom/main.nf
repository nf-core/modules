process VIRUSRECOM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/virusrecom:1.4.0--pyhdfd78af_0':
        'quay.io/biocontainers/virusrecom:1.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(mapping), path(iwic)
    val(is_aligned)
    val(alignment_tool)

    output:
    tuple val(meta), path("${meta.id}/"), emit: results
    tuple val("${task.process}"), val('virusrecom'), eval("virusrecom -h 2>&1 | sed -n 's/.*Version: \\([^ ]*\\).*/\\1/p'"), topic: versions, emit: versions_virusrecom

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def use_iwic      = iwic ? true : false
    def input_command = use_iwic   ? "-iwic ${iwic}"
                      : is_aligned ? "-a ${fasta}"
                      :              "-ua ${fasta}"
    def map_command   = use_iwic ? "" : "-map ${mapping}"
    def at_command    = (!use_iwic && !is_aligned && alignment_tool) ? "-at ${alignment_tool}" : ""

    if (!use_iwic && !is_aligned && !alignment_tool) {
        error "alignment_tool must be specified when input is unaligned and no iwic is provided"
    }
    """
    virusrecom \\
        ${input_command} \\
        ${map_command} \\
        ${at_command} \\
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
