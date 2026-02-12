process MMSEQS_LINCLUST {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    tuple val(meta), path(db_input)

    output:
    tuple val(meta), path("${prefix}/"), emit: db_cluster
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${db_input}" == "${prefix}") {
        error("Input and output names of databases are the same, set prefix in module configuration to disambiguate!")
    }

    """
    mkdir -p ${prefix}
    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_INPUT_PATH_NAME=\$(find -L "${db_input}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        linclust \\
        \$DB_INPUT_PATH_NAME \\
        ${prefix}/${prefix} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}.{0..9}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index

    """
}
