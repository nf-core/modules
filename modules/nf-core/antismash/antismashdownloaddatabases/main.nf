process ANTISMASH_ANTISMASHDOWNLOADDATABASES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "nf-core/antismash:8.0.0"

    output:
    path ("antismash_db"), emit: database
    path ("antismash_dir"), emit: antismash_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download-antismash-databases \\
        --database-dir antismash_db \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "download-antismash-databases --database-dir antismash_db ${args}"

    mkdir antismash_dir
    mkdir antismash_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS
    """
}
