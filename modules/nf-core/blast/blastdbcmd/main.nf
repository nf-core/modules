process BLAST_BLASTDBCMD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(db), val(entry)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    blastdbcmd \\
        -db ${db} \\
        -entry ${entry}
        ${args}
        -out ${meta.id}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastdbcmd -version 2>&1 | sed 's/^.*blastdbcmd: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${meta.id}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastdbcmd -version 2>&1 | sed 's/^.*blastdbcmd: //; s/ .*\$//')
    END_VERSIONS
    """
}

process BLAST_BLASTDBCMD_BATCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(db), path(entry_batch)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    blastdbcmd \\
        -db ${db} \\
        -entry_batch ${entry_batch}
        ${args}
        -out ${meta.id}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastdbcmd -version 2>&1 | sed 's/^.*blastdbcmd: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${meta.id}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastdbcmd -version 2>&1 | sed 's/^.*blastdbcmd: //; s/ .*\$//')
    END_VERSIONS
    """
}
