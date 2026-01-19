process HHSUITE_HHSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hhsuite:3.3.0--py311pl5321h9f068be_13':
        'biocontainers/hhsuite:3.3.0--py311pl5321h9f068be_13' }"

    input:
    tuple val(meta) , path(aln)
    tuple val(meta2), path(hh_db)

    output:
    tuple val(meta), path("*.hhr"), emit: hhr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = aln.getExtension() == "gz" ? true : false
    def aln_name = is_compressed ? aln.getBaseName() : aln
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${aln} > ${aln_name}
    fi

    file=\$(find ${hh_db}/ -type f -name '*_cs219.ffdata' | head -n 1)
    base=\$(basename "\$file")
    db_name=\${base%_cs219.ffdata}

    hhsearch \\
        $args \\
        -cpu $task.cpus \\
        -i ${aln_name} \\
        -d ${hh_db}/\$db_name \\
        -o ${prefix}.hhr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhsearch -h | grep 'HHsearch' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hhr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhsearch -h | grep 'HHsearch' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
