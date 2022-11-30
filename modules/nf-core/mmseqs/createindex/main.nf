process MMSEQS_CREATEINDEX {
    tag "$db"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mmseqs2=14.7e284" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0':
        'quay.io/biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0' }"

    input:
    path db

    output:
    path(db)           , emit: db_indexed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mmseqs createindex \\
        ${db}/${db} \\
        tmp1 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${db}.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
