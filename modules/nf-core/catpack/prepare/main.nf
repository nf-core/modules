process CATPACK_PREPARE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_0'
        : 'biocontainers/cat:6.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(db_fasta)
    path names
    path nodes
    path acc2tax

    output:
    tuple val(meta), path("${prefix}/"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CAT_pack prepare \\
        -n ${task.cpus} \\
        --db_fasta ${db_fasta} \\
        --names ${names} \\
        --nodes ${nodes} \\
        --acc2tax ${acc2tax} \\
        --db_dir ${prefix}/ \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """
}
