process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "bioconda::bakta=1.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.6.0--pyhdfd78af_0' :
        'quay.io/biocontainers/bakta:1.6.0--pyhdfd78af_0' }"

    output:
    path "db.tar.gz"        , emit: db_tar_gz
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bakta_db \\
        download \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --help) 2>&1 | sed 's/.*Version: //g;s/ DOI.*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    echo "bakta_db \\
        download \\
        $args"

    touch db.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --help) 2>&1 | sed 's/.*Version: //g;s/ DOI.*//g')
    END_VERSIONS
    """
}
