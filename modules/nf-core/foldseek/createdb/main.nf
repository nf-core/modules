process FOLDSEEK_CREATEDB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0':
        'biocontainers/foldseek:9.427df8a--pl5321hb365157_0' }"

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
    def dbs = task.ext.dbs ?: "${prefix}"
    """
    mkdir -p ${dbs}
    foldseek \\
        createdb \\
        ${pdb} \\
        ${dbs}/${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbs = task.ext.dbs ?: "${prefix}"

    """
    mkdir -p ${dbs}
    touch ${dbs}/${prefix}
    touch ${dbs}/${prefix}_ca
    touch ${dbs}/${prefix}_ca.dbtype
    touch ${dbs}/${prefix}_ca.index
    touch ${dbs}/${prefix}_h
    touch ${dbs}/${prefix}_h.dbtype
    touch ${dbs}/${prefix}_h.index
    touch ${dbs}/${prefix}_ss
    touch ${dbs}/${prefix}_ss.dbtype
    touch ${dbs}/${prefix}_ss.index
    touch ${dbs}/${prefix}.dbtype
    touch ${dbs}/${prefix}.index
    touch ${dbs}/${prefix}.lookup
    touch ${dbs}/${prefix}.source

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}
