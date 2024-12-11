process MASHMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mashmap:3.1.3--h07ea13f_0':
        'biocontainers/mashmap:3.1.3--h07ea13f_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.paf"), emit: paf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mashmap -q ${fasta} \\
    -r ${reference} \\
    -t ${task.cpus} \\
    -o ${prefix}.paf \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashmap: \$(mashmap --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashmap: \$(mashmap --version 2>&1)
    END_VERSIONS
    """
}
