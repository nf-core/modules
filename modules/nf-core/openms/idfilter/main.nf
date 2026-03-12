process OPENMS_IDFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(id_file), path(filter_file)

    output:
    tuple val(meta), path("*.{idXML,consensusXML}"), emit: filtered
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${id_file.getExtension()}"
    if ("$id_file" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    // Optional filtering via filter_file
    def filter_citerion = task.ext.args2 ?: "-whitelist:peptides"
    def filter = filter_file ? "${filter_citerion} ${filter_file}" : ""

    """
    IDFilter -in $id_file \\
        -out ${prefix}.${suffix} \\
        -threads $task.cpus \\
        $filter \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "${id_file.getExtension()}"
    if ("$id_file" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.${suffix}
    """
}
