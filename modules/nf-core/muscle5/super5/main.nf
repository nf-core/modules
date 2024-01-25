process MUSCLE5_SUPER5 {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/muscle:5.1--h9f5acd7_1':
        'biocontainers/muscle:5.1--h9f5acd7_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    prefix = args.contains('-perm all') ? "${prefix}@" : "${prefix}"
    """
    muscle \\
        -super5 ${fasta} \\
        -output ${prefix}.aln \\
        ${args} \\
        -threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muscle: \$(muscle -version | head -n 1 | cut -d ' ' -f 2 | sed 's/.linux64//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muscle: \$(muscle -version | head -n 1 | cut -d ' ' -f 2 | sed 's/.linux64//')
    END_VERSIONS
    """
}
