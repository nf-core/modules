process MMSEQS_CREATEDB {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0':
        'biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${prefix}/") , emit: db
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    mmseqs \\
        createdb \\
        ${sequence} \\
        ${prefix}/seq \\
        $args \\
        --compressed 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}/seq
    touch ${prefix}/seq.dbtype
    touch ${prefix}/seq.index
    touch ${prefix}/seq.lookup
    touch ${prefix}/seq.source
    touch ${prefix}/seq_h
    touch ${prefix}/seq_h.dbtype
    touch ${prefix}/seq_h.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
