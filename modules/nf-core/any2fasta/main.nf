process ANY2FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/any2fasta:0.8.1--hdfd78af_0'
        : 'biocontainers/any2fasta:0.8.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val("${task.process}"), val('any2fasta'), eval("any2fasta -v 2>&1 | head -1 | sed 's/any2fasta //'"), topic: versions, emit: versions_any2fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    any2fasta \\
        ${args} \\
        ${sequence} \\
        > ${prefix}.fasta
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}.fasta
    """
}
