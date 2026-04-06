process GANON_CLASSIFY {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ganon:2.1.0--py310hab1bfa5_1'
        : 'biocontainers/ganon:2.1.0--py310hab1bfa5_1'}"

    input:
    tuple val(meta), path(fastqs)
    path db

    output:
    tuple val(meta), path("*.tre"), emit: tre
    tuple val(meta), path("*.rep"), emit: report
    tuple val(meta), path("*.one"), emit: one, optional: true
    tuple val(meta), path("*.all"), emit: all, optional: true
    tuple val(meta), path("*.unc"), emit: unc, optional: true
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('ganon'), eval("ganon --version 2>1 | sed 's/.*ganon //g'"), emit: versions_ganon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "--single-reads ${fastqs}" : "--paired-reads ${fastqs}"
    """
    dbprefix=\$(find -L . -name '*.*ibf' | sed 's/\\.h\\?ibf\$//')

    ganon \\
        classify \\
        --db-prefix \${dbprefix%%.*ibf} \\
        ${args} \\
        --threads ${task.cpus} \\
        --output-prefix ${prefix} \\
        ${input} \
        2>&1 | tee ${prefix}.log

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tre
    touch ${prefix}.rep
    touch ${prefix}.one
    touch ${prefix}.all
    touch ${prefix}.unc
    touch ${prefix}.log

    """
}
