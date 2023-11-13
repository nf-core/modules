process FOLDSEEK_CREATEDB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek:8.ef4e960--pl5321hb365157_0':
        'biocontainers/foldseek:8.ef4e960--pl5321hb365157_0' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${meta.id}"), emit: db
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    foldseek \\
        createdb \\
        ${pdb} \\
        ${prefix}/${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p targetDB
    touch ${prefix}
    touch ${prefix}_ca
    touch ${prefix}_ca.dbtype
    touch ${prefix}_ca.index
    touch ${prefix}_h
    touch ${prefix}_h.dbtype
    touch ${prefix}_h.index
    touch ${prefix}_ss
    touch ${prefix}_ss.dbtype
    touch ${prefix}_ss.index
    touch ${prefix}.dbtype
    touch ${prefix}.index
    touch ${prefix}.lookup
    touch ${prefix}.source
    mv ${prefix}* targetDB

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}
