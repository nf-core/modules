process SYLPH_SKETCHGENOMES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(fasta, stageAs: 'genomes/')

    output:
    tuple val(meta), path('*.syldb'), emit: syldb
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls -1 genomes/* > genomes.txt

    sylph sketch \\
        -t ${task.cpus} \\
        ${args} \\
        --gl genomes.txt \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V|awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    echo "${args}"
    touch ${prefix}.syldb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V|awk '{print \$2}')
    END_VERSIONS
    """
}
