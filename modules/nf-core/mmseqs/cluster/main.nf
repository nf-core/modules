process MMSEQS_CLUSTER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0':
        'biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0' }"

    input:
    tuple val(meta), path(db_input)

    output:
    tuple val(meta), path("${prefix}/"), emit: cluster
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$db_input" == "${prefix}") error "Input and output names of databases are the same, set prefix in module configuration to disambiguate!"

    """
    mkdir -p ${prefix}
    DB_INPUT_PATH_NAME=\$(find -L "$db_input/" -name "*.dbtype" | sed 's/\\.dbtype\$//' | head -n 1 )

    mmseqs \\
        cluster \\
        \$DB_INPUT_PATH_NAME \\
        ${prefix}/clust \\
        tmp1 \\
        --threads ${task.cpus} \\
        --compressed 1 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}/clust.{0..9}
    touch ${prefix}/clust.dbtype
    touch ${prefix}/clust.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: /')
    END_VERSIONS
    """
}
