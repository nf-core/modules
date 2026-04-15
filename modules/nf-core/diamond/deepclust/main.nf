process DIAMOND_DEEPCLUST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.24--hf93d47f_0'
        : 'biocontainers/diamond:2.1.24--hf93d47f_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: clusters
    tuple val("${task.process}"), val('diamond'), eval("diamond --version | sed 's/diamond version //g'"), emit: versions_diamond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    diamond \\
        deepclust \\
        ${args} \\
        -p ${task.cpus} \\
        --db ${fasta} \\
        --out ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    touch ${prefix}.tsv
    """
}
