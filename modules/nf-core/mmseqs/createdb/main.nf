process MMSEQS_CREATEDB {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${prefix}/"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = sequence.getExtension() == "gz" ? true : false
    def sequence_name = is_compressed ? sequence.getBaseName() : sequence
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${sequence} > ${sequence_name}
    fi

    mkdir -p ${prefix}

    mmseqs \\
        createdb \\
        ${sequence_name} \\
        ${prefix}/${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index
    touch ${prefix}/${prefix}.lookup
    touch ${prefix}/${prefix}.source
    touch ${prefix}/${prefix}_h
    touch ${prefix}/${prefix}_h.dbtype
    touch ${prefix}/${prefix}_h.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
