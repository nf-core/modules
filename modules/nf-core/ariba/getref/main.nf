process ARIBA_GETREF {
    tag "$db_name"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ariba:2.14.6--py39h67e14b5_3':
        'biocontainers/ariba:2.14.6--py39h67e14b5_3' }"

    input:
    tuple val(meta), val(db_name)

    output:
    tuple val(meta), path("${db_name}.tar.gz"), emit: db
    tuple val("${task.process}"), val('ariba'), eval('ariba version | sed -nE "s/ARIBA version: //p"'), emit: versions_ariba, topic: versions

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
    """

    stub:
    """
    echo "" | gzip > ${db_name}.tar.gz
    """
}
