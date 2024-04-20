process BLAST_UPDATEBLASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta)

    output:
    tuple val(meta), path("${meta.id}"), emit: db
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}; cd ${prefix}

    update_blastdb.pl \\
        ${prefix} \\
        --decompress \\
        ${args}

    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        updateblastdb: \$(update_blastdb.pl -version 2>&1 | tail -n1 | sed 's/[^0-9]*//;')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        updateblastdb: \$(update_blastdb.pl -version 2>&1 | tail -n1 | sed 's/[^0-9]*//;')
    END_VERSIONS
    """
}
