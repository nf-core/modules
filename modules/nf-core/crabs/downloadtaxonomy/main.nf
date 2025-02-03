process CRABS_DOWNLOADTAXONOMY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.0.7--pyhdfd78af_0':
        'biocontainers/crabs:1.0.7--pyhdfd78af_0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("*/nucl_gb.accession2taxid"), emit: accession2taxid
    tuple val(meta), path("*/names.dmp")              , emit: names
    tuple val(meta), path("*/nodes.dmp")              , emit: nodes
    tuple val(meta), path("*")                        , emit: downloadtaxonomy_dir
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    crabs --download-taxonomy \\
        --output ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs --help | grep 'CRABS |' | sed 's/.*CRABS | \\(v[0-9.]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/nucl_gb.accession2taxid
    touch ${prefix}/names.dmp
    touch ${prefix}/nodes.dmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs --help | grep 'CRABS |' | sed 's/.*CRABS | \\(v[0-9.]*\\).*/\\1/')
    END_VERSIONS
    """
}
