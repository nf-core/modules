process GENMAP_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::genmap=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmap:1.3.0--h1b792b2_1' :
        'biocontainers/genmap:1.3.0--h1b792b2_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}") , emit: index
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "$meta.id"

    """
    genmap \\
        index \\
        --fasta-file ${fasta} \\
        --index ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version 2>&1 | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "$meta.id"

    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version 2>&1 | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """
}
