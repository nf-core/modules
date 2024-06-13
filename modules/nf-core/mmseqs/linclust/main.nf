process MMSEQS_LINCLUST {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    tuple val(meta), path(db_input)

    output:
    tuple val(meta), path("${prefix}/"), emit: db_cluster
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$db_input" == "${prefix}") error "Input and output names of databases are the same, set prefix in module configuration to disambiguate!"

    """
    mkdir -p ${prefix}
    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_INPUT_PATH_NAME=\$(find -L "$db_input/" -maxdepth 1 -name "$args2" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        linclust \\
        \$DB_INPUT_PATH_NAME \\
        ${prefix}/${prefix} \\
        tmp1 \\
        $args \\
        --threads ${task.cpus} \\
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
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}.{0..9}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
