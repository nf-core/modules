process MEGAN_DAA2INFO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megan:6.25.9--h9ee0642_0':
        'quay.io/biocontainers/megan:6.25.9--h9ee0642_0' }"

    input:
    tuple val(meta), path(daa)
    val(megan_summary)

    output:
    tuple val(meta), path("*.txt.gz"), emit: txt_gz
    tuple val(meta), path("*.megan") , emit: megan, optional: true
    tuple val("${task.process}"), val('megan'), eval("daa2info 2>&1 | sed '/version/!d;s/.*version //;s/, .*//'"), emit: versions_megan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def summary = megan_summary ? "-es ${prefix}.megan" : ""

    """
    daa2info \\
        -i ${daa} \\
        -o ${prefix}.txt.gz \\
        ${summary} \\
        ${args}
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def megan_cmd = megan_summary ? "touch ${prefix}.megan" : ""

    """
    echo "" | gzip > ${prefix}.txt.gz
    ${megan_cmd}
    """
}
