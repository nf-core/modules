process BLAST_UPDATEBLASTDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.17.0--h66d330f_0':
        'biocontainers/blast:2.17.0--h66d330f_0' }"

    input:
    tuple val(meta), val(name)

    output:
    tuple val(meta), path(prefix), emit: db
    tuple val("${task.process}"), val("updateblastdb"), eval("update_blastdb.pl -version 2>&1 | tail -n1 | rev | cut -f1 -d ' ' | rev"), topic: versions, emit: versions_updateblastdb

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}; cd ${prefix}

    update_blastdb.pl \\
        ${name} \\
        ${args}

    cd ..

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${name}.ndb

    """
}
