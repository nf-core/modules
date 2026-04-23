process KMCP_PROFILE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmcp:0.9.4--h9ee0642_0'
        : 'quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0'}"

    input:
    tuple val(meta), path(search_results)
    path db
    val level

    output:
    tuple val(meta), path("*.profile"), emit: profile
    tuple val("${task.process}"), val('kmcp'), eval("kmcp version 2>&1 | sed 's/^.*kmcp v//'"), emit: versions_kmcp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxid = level == 'species' ? "-T \$(find -L . -name '*map' -type f)" : ""
    def taxdmp = level == 'species' ? "-X \$(find -L . -name 'nodes.dmp' -type f | sed 's#nodes.dmp\$##g')" : ""
    if (!['species', 'strain', 'assembly'].contains(level)) {
        error("[KMCP_PROFILE] ERROR: --level must be one of 'species', 'strain' or 'assembly'")
    }
    """
    ## Input validation checks, as require files come via a directory in an input channel (can't come through separate files)
    if [[ "${level}" == "species" && \$(find -L . -name '*map' -type f) == "" ]]; then
        echo "[KMCP_PROFILE] ERROR: --level 'species' requires a reference to taxid map file to be provided in the database directory. These files are not detected."
    fi

    if [[ "${level}" == "species" && \$(find -L . -name 'nodes.dmp' -type f | sed 's#nodes.dmp\$##g') == "" ]]; then
        echo "[KMCP_PROFILE] ERROR: --level 'species' requires taxdump files to be provided in the database directory. These files are not detected."
    fi

    kmcp \\
        profile \\
        --level ${level} \\
        -j ${task.cpus} \\
        -o ${prefix}.profile \\
        ${args} \\
        ${taxdmp} \\
        ${taxid} \\
        ${search_results}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!['species', 'strain', 'assembly'].contains(level)) {
        error("[KMCP_PROFILE] ERROR: --level must be one of 'species', 'strain' or 'assembly'")
    }
    """
    touch ${prefix}.profile
    """
}
