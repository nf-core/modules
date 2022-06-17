process ARIBA_GETREF {
    tag "$db_name"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ariba=2.14.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ariba:2.14.6--py39h67e14b5_3':
        'quay.io/biocontainers/ariba:2.14.6--py39h67e14b5_3' }"

    input:
    val(db_name)

    output:
    tuple path("${db_name}.tar.gz"), emit: db
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Download, format database, and tarball it
    ariba \\
        getref \\
        ${db_name} \\
        ${db_name}

    ariba \\
        prepareref \\
        -f ${db_name}.fa \\
        -m ${db_name}.tsv \\
        ${db_name}

    tar -zcvf ${db_name}.tar.gz ${db_name}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
