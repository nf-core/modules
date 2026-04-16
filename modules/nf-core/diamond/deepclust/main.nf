process DIAMOND_DEEPCLUST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.24--hf93d47f_0'
        : 'biocontainers/diamond:2.1.24--hf93d47f_0'}"

    input:
    tuple val(meta), path(fasta)
    val save_aln

    output:
    tuple val(meta), path("*.tsv"), emit: clusters
    tuple val(meta), path("*.aln"), optional: true, emit: alignment
    tuple val("${task.process}"), val('diamond'), eval("diamond --version | sed 's/.* //g'"), emit: versions_diamond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def aln_out = save_aln ? "--aln-out ${prefix}.aln" : ''
    """
    diamond \\
        deepclust \\
        ${args} \\
        --memory-limit  "${task.memory.toGiga()}G" \\
        -p ${task.cpus} \\
        --db ${fasta} \\
        --out ${prefix}.tsv \\
        ${aln_out}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def aln_out = save_aln ? "touch ${prefix}.aln" : ''
    """
    echo "${args}"
    touch ${prefix}.tsv
    ${aln_out}
    """
}
